#! /usr/bin/env python

##
# Python runner for orbit visualization pipeline
# @author Kit Kennedy
#

import time
import os.path
import json
from datetime import datetime, timedelta
import sys
import argparse
import numpy as np
from copy import deepcopy
from collections import OrderedDict
import multiprocessing as mp


#  local repo includes. todo:  make this less hackey

# add this main tester so that if we're running the gp remotely (e.g. in circinus sim), we don't try to reinclude a bunch of modules that are already present at the remote
if __name__ == "__main__":
    sys.path.append ('..')
from circinus_tools  import time_tools as tt
from circinus_tools  import io_tools
from circinus_tools.scheduling.routing_objects import DataMultiRoute
from circinus_tools.scheduling.io_processing import SchedIOProcessor
from gp_tools.gp_plotting import GPPlotting
import gp_tools.gp_route_selection_v1 as gprsv1
import gp_tools.gp_route_selection_v2 as gprsv2
from gp_tools.gp_activity_scheduling_separate import GPActivitySchedulingSeparate
from gp_tools.gp_activity_scheduling_coupled import GPActivitySchedulingCoupled
from gp_tools.gp_metrics import GPMetrics
import gp_tools.gp_rs_execution as gp_rs_exec
import gp_tools.gp_general_tools as gp_gen

import runner_gp_helper_pickle as pickle_helper 
import runner_gp_helper_other as other_helper
import runner_gp_helper_output as output_helper

from circinus_tools import debug_tools

REPO_BASE = os.path.abspath(os.pardir)  # os.pardir aka '..'

OUTPUT_JSON_VER = '0.2'

def print_verbose(string,verbose=False):
    if verbose:
        print(string)


class GlobalPlannerRunner:
    """easy interface for running the global planner scheduling algorithm"""

    def __init__(self, gp_params):
        """initializes based on parameters

        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """

        self.params = gp_params
        self.scenario_params = self.params['orbit_prop_params']['scenario_params']
        self.sat_params = self.params['orbit_prop_params']['sat_params']
        self.obs_params = self.params['orbit_prop_params']['obs_params']
        self.sat_orbit_params = self.params['orbit_prop_params']['sat_orbit_params']
        self.obs_params = self.params['orbit_prop_params']['obs_params']
        self.pickle_params = self.params['gp_general_params']['pickle_params']
        self.other_params = self.params['gp_other_params']
        self.plot_params = self.params['gp_general_params']['plot_params']
        self.general_other_params = self.params['gp_general_params']['other_params']
        self.rs_general_params = gp_params['gp_general_params']['route_selection_general_params']
        self.rs_v2_params = gp_params['gp_general_params']['route_selection_params_v2']
        self.as_params = self.params['gp_general_params']['activity_scheduling_params']
        self.gp_agent_ID = gp_params['gp_instance_params']['gp_agent_ID']
        self.gp_inst_planning_params = self.params['gp_instance_params']['planning_params']
        self.gp_inst_as_params = self.params['gp_instance_params']['activity_scheduling_params']
        self.io_proc =SchedIOProcessor(self.params)
        self.gp_plot =GPPlotting( self.params)

        # this data route ID should be unique across all data routes (and the datamultiroutes encasing them).  this is important because data route IDs are used for indexing in many places, not just in the global planner.
        self.initial_dr_uid = gp_params['gp_instance_params']['initial_gp_route_indx']

        # it's no good if the planning window (the activities to select) goes outside of the scenario window (the possible activities, eclipse windows)
        if self.gp_inst_planning_params['planning_start_dt'] < self.scenario_params['start_utc_dt']:
            raise RuntimeError("GP instance start time (%s) is less than scenario start time (%s)"%(self.gp_inst_planning_params['planning_start_dt'],self.scenario_params['start_utc']))
        if self.gp_inst_planning_params['planning_end_obs_dt'] > self.scenario_params['end_utc_dt']:
            raise RuntimeError("GP instance obs,xlnk end time (%s) is greater than scenario end time (%s)"%(self.gp_inst_planning_params['planning_end_obs_dt'],self.scenario_params['end_utc']))
        if self.gp_inst_planning_params['planning_end_xlnk_dt'] > self.scenario_params['end_utc_dt']:
            raise RuntimeError("GP instance obs,xlnk end time (%s) is greater than scenario end time (%s)"%(self.gp_inst_planning_params['planning_end_xlnk_dt'],self.scenario_params['end_utc']))
        if self.gp_inst_planning_params['planning_end_dlnk_dt'] > self.scenario_params['end_utc_dt']:
            raise RuntimeError("GP instance dlnk end time (%s) is greater than scenario end time (%s)"%(self.gp_inst_planning_params['planning_end_dlnk_dt'],self.scenario_params['end_utc']))

    



    def run_activity_scheduling( self, routes_by_obs,existing_route_data,ecl_winds,verbose=False):
        gp_as = GPActivitySchedulingSeparate ( self.params)

        # flatten the list of all routes, which currently has nested lists for each observation
        routes_flat = [rt for rts in routes_by_obs.values() for rt in rts]

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        existing_routes = existing_route_data.get('existing_routes',[])
        utilization_by_existing_route_id = existing_route_data.get('utilization_by_existing_route_id',{})

        print_verbose('make activity scheduling model',verbose)
        gp_as.make_model (routes_flat, existing_routes, utilization_by_existing_route_id, ecl_winds,verbose = verbose)
        stats =gp_as.get_stats (verbose = verbose)
        print_verbose('solve activity scheduling',verbose)
        t_a = time.time()
        gp_as.solve ()
        t_b = time.time()
        # gp_as.print_sol ()
        print_verbose('extract_routes',verbose)
        #  make a copy of the windows in the extracted routes so we don't mess with the original objects ( just to be extra careful)
        scheduled_routes,all_updated_routes = gp_as.extract_utilized_routes ( verbose  = False)
        energy_usage,data_usage = gp_as.extract_resource_usage(  decimation_factor =1)
        # look at reasons routes weren't scheduled
        extract_reasoning = False
        if extract_reasoning:
            gp_as.extract_schedule_reasoning( routes_flat, verbose = verbose)

        time_elapsed = t_b-t_a

        return  scheduled_routes, all_updated_routes, energy_usage, data_usage

    def run_activity_scheduling_coupled( self, obs_winds,dlnk_winds_flat,xlnk_winds,ecl_winds,verbose=False):
        gp_as = GPActivitySchedulingCoupled( self.params)

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        print_verbose('make activity scheduling (coupled) model',verbose)
        gp_as.make_model ( obs_winds,dlnk_winds_flat,xlnk_winds,ecl_winds,verbose = verbose)
        stats =gp_as.get_stats (verbose = verbose)
        print_verbose('solve activity scheduling (coupled)',verbose)
        t_a = time.time()
        gp_as.solve ()
        t_b = time.time()
        # gp_as.print_sol ()

        print_verbose('extract_routes',verbose)
        #  make a copy of the windows in the extracted routes so we don't mess with the original objects ( just to be extra careful)
        routes = gp_as.extract_utilized_routes ( copy_routes = True, verbose  = False)
        energy_usage,data_usage = gp_as.extract_resource_usage(  decimation_factor =1)

        # for dmr in routes:
        #     for dr in dmr.data_routes:
        #         print_verbose(dr,verbose)

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        time_elapsed = t_b-t_a

        return  routes, energy_usage, data_usage


    def run_route_selection_v2_step1( self,obs_winds,dlnk_winds_flat,xlnk_winds,existing_route_data,latest_dr_uid,verbose=False):
        # todo: it would be nice at some point to incorporate the use of existing routes into here.  but that's definitely for future work.

        print_verbose ('nominal route selection v2 step 1',verbose)

        obs_indx =0
        # dict of all routes, with obs as key
        routes_by_obs ={}
        indices_by_obs ={}
        all_stats =[]
        route_times_s =[]

        wind_ids_seen = set()
        def filter_obs(obs_winds):
            obs_winds_filt = [[] for sat_indx in range( self.sat_params['num_sats'])]

            for sat_indx in range( self.sat_params['num_sats']):
                for sat_obs_index,obs_wind in enumerate(obs_winds[sat_indx]):
                    #  filter out any redundant windows created by a possible earlier copying somewhere
                    if obs_wind.window_ID in wind_ids_seen:
                        continue

                    # filter for observations that come between end of fixed planning time and end of planning window.  the fixed planning time is the time up to which we are saying no new activity windows are allowed to be scheduled.  before that only existing routes may have activity window scheduled.
                    # if obs_wind.start >= self.gp_inst_planning_params['planning_fixed_end_dt'] and obs_wind.end <= self.gp_inst_planning_params['planning_end_obs_dt']:
                    if gp_gen.wind_in_planning_window(self,obs_wind):
                        obs_winds_filt[sat_indx].append(obs_wind)
                        wind_ids_seen.add(obs_wind)
            return obs_winds_filt


        obs_winds_filt = filter_obs(obs_winds)


        # this may seem mildy pointless at first, but what we're doing is making a flat lookup table from the hash of wind to the wind object itself. We'll need this to restore the original winds later after running in parallel.
        obs_winds_dict = {wind:wind for sat_indx in range (self.sat_params['num_sats']) for wind in obs_winds_filt[sat_indx] }
        dlnk_winds_dict = {wind:wind for sat_indx in range (self.sat_params['num_sats']) for wind in dlnk_winds_flat[sat_indx] }
        xlnk_winds_dict = {wind:wind for sat_indx in range (self.sat_params['num_sats']) for xsat_indx in range (self.sat_params['num_sats']) for wind in xlnk_winds[sat_indx][xsat_indx]}

        num_obs = sum(len(obs_winds_filt[sat_indx]) for sat_indx in range (self.sat_params['num_sats']))

        t_a = time.time()

        def org_outputs(obs_output,obs_indx,sat_indx,sat_obs_index):
            #  unpack items from the output tuple
            obs_wind,routes,time_elapsed = obs_output
            routes_by_obs[obs_wind] = routes
            # todo: should probably make use of this...
            all_stats.append ( None)
            route_times_s.append ( time_elapsed)
            indices_by_obs[obs_wind] = [obs_indx,sat_indx,sat_obs_index ]


        #  run in parallel mode
        if self.rs_general_params['run_rs_parallel']:
            #  make the inputs
            all_obs_inputs = []
            obs_indx = 0
            for sat_indx in range( self.sat_params['num_sats']):
                for sat_obs_index,obs_wind in  enumerate(obs_winds_filt[sat_indx]):
                    all_obs_inputs.append((obs_wind,obs_indx,sat_indx,sat_obs_index))
                    obs_indx += 1

            p = mp.Pool(self.rs_general_params['num_parallel_workers'])
            obs_outputs_list = p.map(gp_rs_exec.ParallelRSWorkerWrapper(gprsv2.GPDataRouteSelection,self.params,dlnk_winds_flat,xlnk_winds,num_obs, verbose=verbose), all_obs_inputs)
            p.terminate()

            #  unpack all of the outputs from the parallel runs list
            for output_indx,obs_output in  enumerate(obs_outputs_list):
                #  same index into the inputs as the outputs
                _,obs_indx,sat_indx,sat_obs_index = all_obs_inputs[output_indx]
                org_outputs(obs_output,obs_indx,sat_indx,sat_obs_index)

            # because we duplicated dlnk_winds_flat and xlnk_winds for each parallel worker, need to restore the original window objects in all the routes
            gp_rs_exec.restore_original_wind_obj(routes_by_obs,obs_winds_dict,dlnk_winds_dict,xlnk_winds_dict)

        #  run in serial mode
        else:
            gp_rs = gprsv2.GPDataRouteSelection ( self.params)

            for sat_indx in range( self.sat_params['num_sats']):
                for  sat_obs_index, obs_wind in  enumerate ( obs_winds_filt[sat_indx]):
                    if verbose:
                        gp_rs_exec.print_pre_run_msg(obs_wind, obs_indx, sat_indx, sat_obs_index,num_obs)

                    obs_output = gp_rs_exec.run_rs_step1(gp_rs,obs_wind,dlnk_winds_flat,xlnk_winds,verbose = False)

                    org_outputs(obs_output,obs_indx,sat_indx,sat_obs_index)

                    obs_indx +=1

                    # print_verbose ('len(routes)',verbose)
                    # print_verbose (len(routes),verbose)
                    # print_verbose ('obs.data_vol, total dlnk dv, ratio',verbose)
                    # print_verbose (obs.data_vol,sum(dr.data_vol for dr in routes),sum(dr.data_vol for dr in routes)/obs.data_vol,verbose)
                    # print_verbose ('min latency, ave latency, max latency',verbose)
                    # latencies = [(rt.route[-1].end - rt.route[0].end).total_seconds()/60 for rt in  routes]
                    # if len(latencies) > 0:
                    #     print_verbose (np.min( latencies),np.mean( latencies),np.max( latencies),verbose)
                    # else:
                    #     print_verbose('no routes found',verbose)

        t_b = time.time()
        time_elapsed = t_b-t_a


        # explicitly validate routes
        for routes in routes_by_obs.values():
            for dr in routes:
                #  set the data route ID now, because they were not set uniquely if we were running in parallel
                dr.set_id(self.gp_agent_ID,latest_dr_uid)
                latest_dr_uid += 1

                # try:
                dr.validate()
                # except:
                #     raise Exception("Couldn't handle dr %s for obs %s (indices %s)"%(dr,dr.get_obs(),indices_by_obs[dr.get_obs()]))

        if verbose:
            print_verbose('total_rs_step1_runtime',verbose)
            print_verbose(time_elapsed,verbose)
            valid = len(route_times_s) > 0
            if valid:
                print_verbose('runtime per rs iteration',verbose)
                print_verbose('ave: %fs'%(np.mean(route_times_s)),verbose)
                print_verbose('std: %fs'%(np.std(route_times_s)),verbose)

        return routes_by_obs,all_stats,route_times_s,obs_indx, latest_dr_uid



    def run_route_selection_v2_step2(self,routes_by_obs,existing_route_data,verbose=False):

        routes_by_obs_filt = gp_gen.filt_routes_by_obs(self,routes_by_obs)

        existing_routes = existing_route_data.get('existing_routes',[])
        utilization_by_existing_route_id = existing_route_data.get('utilization_by_existing_route_id',{})

        print_verbose('num routes',verbose)
        print_verbose(sum(len(rts) for rts in routes_by_obs_filt.values()),verbose)
        print_verbose('num routes existing',verbose)
        print_verbose(len(existing_routes),verbose)

        #  go ahead and throw in the existing routes as well.  note that in the general case existing routes might be outside of the filter window.  this is okay because the algorithms below should not depend on data routes being within the filter times.  if route selection step two happens to not choose one or more of the existing routes, that's okay - we'll add them back in at activity scheduling.  however we do want to feed them into step two so that we don't choose new routes that are basically the equivalent of the existing routes just because we don't know that the existing routes exist.
        for rt in existing_routes:
            obs = rt.get_obs()
            # note this comparison is based on window ID,  so if the obs window object is duplicated it should still work okay
            routes_by_obs_filt.setdefault(obs,[]).append(rt)
        existing_routes_set = set(existing_routes)

        gp_met = GPMetrics(self.params)
        print_verbose('Assess route overlap pre RS step 2',verbose)
        overlap_cnt_by_route,stats_rs2_pre = gp_met.assess_route_overlap( routes_by_obs_filt,verbose=verbose)
        # print_verbose(time_elapsed,verbose)

        gp_rs = gprsv2.GPDataRouteSelection ( self.params)

        t_a = time.time()
        selected_rts_by_obs = gp_rs.run_step2(routes_by_obs_filt,overlap_cnt_by_route,existing_routes_set)
        t_b = time.time()
        time_elapsed = t_b-t_a
        print_verbose('RS step2 time_elapsed %f'%(time_elapsed),verbose)

        print_verbose('Assess route overlap post RS step 2',verbose)
        overlap_cnt_by_route,stats_rs2_post = gp_met.assess_route_overlap( selected_rts_by_obs,verbose=verbose)


        for rts in selected_rts_by_obs.values():
            for dmr in rts:
                dmr.validate()

        return selected_rts_by_obs,stats_rs2_pre,stats_rs2_post

    def run_route_selection(self,obs_winds,dlnk_winds_flat,xlnk_winds,ecl_winds,existing_route_data,window_uid,latest_dr_uid,verbose=False):


        #################################
        #  route selection step 1
        #################################

        sel_routes_by_obs = None
        pas_a = None

        # If we need output from step 1
        run_step_1 = not self.other_params['rs_s2_pickle_input'] and not self.other_params['as_pickle_input']
        if run_step_1:

            #  if  we are loading from file, do that
            if self.other_params['rs_s1_pickle_input']:
                print_verbose('Unpickling route selection step one stuff',verbose)
                routes_by_obs,all_stats,route_times_s,obs_indx,ecl_winds,window_uid = pickle_helper.unpickle_rtsel_s1_stuff(self)
                pas_a = time.time()

            #  otherwise run route selection step 1
            else:
                print_verbose('Run route selection step 1',verbose)
                routes_by_obs,all_stats,route_times_s, obs_indx, latest_dr_uid  =  self.run_route_selection_v2_step1(obs_winds,dlnk_winds_flat,xlnk_winds,existing_route_data,latest_dr_uid,verbose=self.rs_v2_params['verbose_step1'])
                # routes_by_obs,all_stats,route_times_s, obs_indx, weights_tups  =  self.run_test_route_selection(obs_winds,dlnk_winds_flat,xlnk_winds)

            #  pickle before step 2 because step 2 doesn't take that long
            if self.pickle_params['pickle_route_selection_step1_results']:
                pickle_helper.pickle_rtsel_s1_stuff(self,routes_by_obs,all_stats,route_times_s,obs_indx,ecl_winds,window_uid)
        else:
            print_verbose('Skipping route selection step one stuff',verbose)


        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        #################################
        #  route selection step 2
        #################################


        # If we need output from step 2
        run_step_2 = not self.other_params['as_pickle_input']
        if run_step_2:
            if self.other_params['rs_s2_pickle_input']:
                print_verbose('Unpickling route selection step two stuff',verbose)
                sel_routes_by_obs,ecl_winds,window_uid,stats_rs2_pre,stats_rs2_post,latest_dr_uid = pickle_helper.unpickle_rtsel_s2_stuff(self)
                pas_a = time.time()
            else:
                rs_s1_found_obs= len(routes_by_obs.keys()) > 0
                if rs_s1_found_obs:
                    print_verbose('num_obs_calced',verbose)
                    print_verbose(len(route_times_s),verbose)
                    print_verbose('np.mean(route_times_s)',verbose)
                    print_verbose(np.mean(route_times_s),verbose)
                    print_verbose('np.max(route_times_s)',verbose)
                    print_verbose(np.max(route_times_s),verbose)
                    print_verbose('np.std(route_times_s)',verbose)
                    print_verbose(np.std(route_times_s),verbose)

                passthru = False
                # passthru is what you use if you just want to feed ALL the data routes from S1 to activity scheduling
                if passthru:
                    stats_rs2_pre = stats_rs2_post = []
                    routes_by_obs_filt = gp_gen.filt_routes_by_obs(self,routes_by_obs)  # may not have been filtered already, if we're loading from a pickle
                    sel_routes_by_obs = {obs:[DataMultiRoute(dr.ID,data_routes=[dr]) for dr in rts] for obs,rts in 
                    routes_by_obs_filt.items()}
                else:
                    # run step 2. todo:  move this elsewhere
                    sel_routes_by_obs,stats_rs2_pre,stats_rs2_post = self.run_route_selection_v2_step2(routes_by_obs,existing_route_data,verbose)

            if self.pickle_params['pickle_route_selection_step2_results']:
                pickle_helper.pickle_rtsel_s2_stuff(self,sel_routes_by_obs,ecl_winds,window_uid,stats_rs2_pre,stats_rs2_post,latest_dr_uid)
        else:
            print_verbose('Skipping route selection step two stuff',verbose)

        #################################
        # route selection output stage
        #################################

        print_verbose('route selection output stage',verbose)

        if self.rs_general_params['plot_route_selection_results']:
            # todo: currently broken - update
            output_helper.plot_route_selection_results (self,sel_routes_by_obs,dlnk_winds_flat,xlnk_winds_flat,num_obs_to_plot = 5)

        return sel_routes_by_obs,ecl_winds,latest_dr_uid,window_uid,pas_a

    def run( self, existing_route_data, verbose=False):

        print_verbose('planning_start_dt: %s'%(self.gp_inst_planning_params['planning_start_dt']),verbose)
        print_verbose('planning_end_obs_dt: %s'%(self.gp_inst_planning_params['planning_end_obs_dt']),verbose)
        print_verbose('planning_end_xlnk_dt: %s'%(self.gp_inst_planning_params['planning_end_xlnk_dt']),verbose)
        print_verbose('planning_end_dlnk_dt: %s'%(self.gp_inst_planning_params['planning_end_dlnk_dt']),verbose)

        #################################
        #  parse inputs, if desired
        #################################

        load_windows = (not self.other_params['rs_s1_pickle_input'] and not self.other_params['rs_s2_pickle_input']) or self.gp_inst_as_params['plot_activity_scheduling_results'] or self.rs_general_params['plot_route_selection_results']

        if load_windows:
            print_verbose('Load files',verbose)

            # parse the inputs into activity windows
            window_uid = 0
            print_verbose('Load obs',verbose)
            obs_winds, window_uid =self.io_proc.import_obs_winds(window_uid)
            print_verbose('Load dlnks',verbose)
            dlnk_winds, dlnk_winds_flat, window_uid =self.io_proc.import_dlnk_winds(window_uid)
            print_verbose('Load xlnks',verbose)
            xlnk_winds, xlnk_winds_flat, window_uid =self.io_proc.import_xlnk_winds(window_uid)
            # note: this import is currently done independently from circinus constellation sim. If we ever need to share knowledge about ecl winds between the two, will need to make ecl winds an input from const sim
            print_verbose('Load ecl',verbose)
            ecl_winds, window_uid =self.io_proc.import_eclipse_winds(window_uid)

            # with open('temp.pkl','wb') as f:
            #     pickle.dump( {'params': self.params},f)

            # important to check this because window unique IDs are used as hashes in dictionaries in the scheduling code
            print_verbose('Validate windows',verbose)
            other_helper.validate_unique_windows(self,obs_winds,dlnk_winds_flat,xlnk_winds,ecl_winds)

            # todo:  probably ought to delete the input times and rates matrices to free up space

            print_verbose('In windows loaded from file:',verbose)
            print_verbose('obs_winds',verbose)
            print_verbose(sum([len(p) for p in obs_winds]),verbose)
            print_verbose('dlnk_win',verbose)
            print_verbose(sum([len(p) for p in dlnk_winds]),verbose)
            print_verbose('xlnk_win',verbose)
            print_verbose(sum([len(xlnk_winds[i][j]) for i in  range( self.sat_params['num_sats']) for j in  range( self.sat_params['num_sats']) ]),verbose)
        else:
            obs_winds = []
            dlnk_winds_flat = []
            xlnk_winds = []
            ecl_winds = []
            window_uid = 0

        #################################
        #  route selection stage
        #################################

        latest_dr_uid = self.initial_dr_uid

        # pas = planning and scheduling
        pas_a = time.time()

        run_rs = not self.as_params['run_coupled_rs_as']

        pas_a_new = None
        if run_rs:
            sel_routes_by_obs,ecl_winds,latest_dr_uid,window_uid,pas_a_new = self.run_route_selection(obs_winds,dlnk_winds_flat,xlnk_winds,ecl_winds,existing_route_data,window_uid,latest_dr_uid,verbose)

        # debug_tools.debug_breakpt()


        if pas_a_new:
            pas_a = pas_a_new

        #################################
        #  activity scheduling stage
        #################################

        if not self.as_params['run_activity_scheduling']:
            return [],[]

        if self.other_params['as_pickle_input']:
            sel_routes_by_obs,ecl_winds,scheduled_routes,energy_usage,data_usage, window_uid,latest_dr_uid = pickle_helper.unpickle_actsc_stuff(self)
        else:
            run_coupled_rs_as = self.as_params['run_coupled_rs_as']

            # make a func so it only gets called if needed
            def found_routes(): 
                return any([len(rts) >0 for rts in sel_routes_by_obs.values()])

            if not run_coupled_rs_as and found_routes():
                #  to protect against the weird case where we didn't find any routes ( shouldn't happen, unless we're at the very end of the simulation, or you're trying to break things)
                # if self.as_params['validate_unique_wind_objects']:  
                #     other_helper.validate_unique_window_objects(self,sel_routes_by_obs,existing_route_data)

                scheduled_routes,all_updated_routes,energy_usage,data_usage = self.run_activity_scheduling(sel_routes_by_obs,existing_route_data,ecl_winds,verbose)

            # run coupled route selection/act sched solver (slow, optimal)
            elif run_coupled_rs_as:
                scheduled_routes,energy_usage,data_usage = self.run_activity_scheduling_coupled(obs_winds,dlnk_winds_flat,xlnk_winds_flat,ecl_winds,verbose)

                #  this is not used at all in the coupled activity scheduling solution, because it currently can't be used in the constellation simulation ( not implemented)
                all_updated_routes = []

                # we didn't run RS, so there are no "selected routes", just scheduled
                sel_routes_by_obs = {}
            else:
                scheduled_routes,all_updated_routes,energy_usage,data_usage = ([],[],None,None)
                print_verbose('No routes were found in route selection; not running activity selection',verbose)

        # if we are saving to file, do that
        if self.pickle_params['pickle_act_scheduling_results']:
            pickle_helper.pickle_actsc_stuff(self,sel_routes_by_obs,ecl_winds,scheduled_routes,energy_usage,data_usage, window_uid,latest_dr_uid)


        pas_b = time.time()
        total_plan_and_sched_runtime = pas_b - pas_a

        #  note that these calculations may be off when running the constellation simulation. don't rely on these for the constellation simulation.  they could be off because there might be duplicate/copy observation windows
        metrics_plot_inputs = output_helper.calc_activity_scheduling_results (self,obs_winds,dlnk_winds_flat,sel_routes_by_obs,scheduled_routes, energy_usage)

        print_verbose('total_plan_and_sched_runtime (warning: may include (un)pickling time and RS plot output)',verbose)
        print_verbose("%.2f seconds"%(total_plan_and_sched_runtime),verbose)

        #################################
        #  Activity scheduling output stage
        #################################

        print_verbose('activity scheduling output stage',verbose)

        if self.gp_inst_as_params['plot_activity_scheduling_results']:
            output_helper.plot_activity_scheduling_results(self,(obs_winds,dlnk_winds_flat,xlnk_winds_flat),sel_routes_by_obs,scheduled_routes,energy_usage,data_usage,ecl_winds,metrics_plot_inputs)


        # if you want to see windows from RS output...
        # sel_routes_flat = [dr for rts in sel_routes_by_obs.values() for dr in rts]
        # (sel_obs_winds_flat, sel_dlnk_winds_flat, sel_xlnk_winds_flat, link_info_by_wind, route_ids_by_wind) = self.io_proc.extract_flat_windows (sel_routes_flat)
        # outputs= self.io_proc.make_sat_history_outputs (sel_obs_winds_flat, sel_xlnk_winds_flat, sel_dlnk_winds_flat, link_info_by_wind)

        (sched_obs_winds_flat, sched_dlnk_winds_flat, sched_xlnk_winds_flat, link_info_by_wind, route_ids_by_wind) = self.io_proc.extract_flat_windows (scheduled_routes)
        viz_outputs= self.io_proc.make_sat_history_outputs (sched_obs_winds_flat, sched_xlnk_winds_flat, sched_dlnk_winds_flat, link_info_by_wind)


        return scheduled_routes,all_updated_routes,viz_outputs,latest_dr_uid


class PipelineRunner:

    def run(self, data,verbose=False):
        """

        """

        # deepcopy here because we're making changes for components that this function calls, and don't want to accidentally somehow step on toes somewhere down the call stack (before this function was called)
        orbit_prop_inputs = deepcopy( data['orbit_prop_inputs'])
        orbit_link_inputs = data['orbit_link_inputs']
        gp_general_params_inputs = data['gp_general_params_inputs']
        gp_instance_params_inputs = data['gp_instance_params_inputs']
        data_rates_inputs = data['data_rates_inputs']
        file_params = data.get ('file_params',{})

        gp_params = {}
        orbit_prop_params = orbit_prop_inputs
        orbit_link_params = orbit_link_inputs
        gp_general_params = gp_general_params_inputs
        gp_instance_params = gp_instance_params_inputs
        data_rates_params = data_rates_inputs
        gp_other_params = {}
        gp_other_params['new_pickle_file_name_pre']  = file_params.get ('new_pickle_file_name_pre' ,'default_pickle')
        gp_other_params['rs_s1_pickle_input']  = data['rs_s1_pickle']
        gp_other_params['rs_s2_pickle_input']  = data['rs_s2_pickle']
        gp_other_params['as_pickle_input']  = data['as_pickle']

        if data['rs_s1_pickle'] and data['rs_s2_pickle']:
            raise Exception('Should only specify 1 input pickle for route selection')

        if orbit_prop_inputs['version'] == "0.5":
            # do some useful transformations while preserving the structure of the inputs ( important for avoiding namespace clashes)
            orbit_prop_inputs['scenario_params']['start_utc_dt'] = tt.iso_string_to_dt ( orbit_prop_inputs['scenario_params']['start_utc'])
            orbit_prop_inputs['scenario_params']['end_utc_dt'] = tt.iso_string_to_dt ( orbit_prop_inputs['scenario_params']['end_utc'])
            orbit_prop_inputs['sat_params']['num_sats'] = orbit_prop_inputs['sat_params']['num_satellites']
            orbit_prop_inputs['gs_params']['num_gs'] = orbit_prop_inputs['gs_params']['num_stations']
            orbit_prop_inputs['sat_params']['pl_data_rate'] = orbit_prop_inputs['sat_params']['payload_data_rate_Mbps']
            # orbit_prop_inputs['sat_orbit_params'], dummy = io_tools.unpack_sat_entry_list( orbit_prop_inputs['sat_orbit_params'],force_duplicate =  True)
            orbit_prop_inputs['sat_params']['power_params'], all_sat_ids1 = io_tools.unpack_sat_entry_list( orbit_prop_inputs['sat_params']['power_params'])
            orbit_prop_inputs['sat_params']['data_storage_params'], all_sat_ids2 = io_tools.unpack_sat_entry_list( orbit_prop_inputs['sat_params']['data_storage_params'])
            orbit_prop_inputs['sat_params']['initial_state'], all_sat_ids3 = io_tools.unpack_sat_entry_list( orbit_prop_inputs['sat_params']['initial_state'])

            #  check if  we saw the same list of satellite IDs from each unpacking. if not that's a red flag that the inputs could be wrongly specified
            if all_sat_ids1 != all_sat_ids2 or all_sat_ids1 != all_sat_ids3:
                raise Exception('Saw differing sat ID orders')

            #  grab the list for satellite ID order.  if it's "default", we will create it and save it for future use here
            sat_id_order=orbit_prop_inputs['sat_params']['sat_id_order']
            #  make the satellite ID order. if the input ID order is default, then will assume that the order is the same as all of the IDs found in the power parameters
            sat_id_order = io_tools.make_and_validate_sat_id_order(sat_id_order,orbit_prop_inputs['sat_params']['num_sats'],all_sat_ids1)
            orbit_prop_inputs['sat_params']['sat_id_order'] = sat_id_order

            orbit_prop_inputs['sat_params']['power_params_sorted'] = io_tools.sort_input_params_by_sat_IDs(orbit_prop_inputs['sat_params']['power_params'],sat_id_order)
            orbit_prop_inputs['sat_params']['data_storage_params_sorted'] = io_tools.sort_input_params_by_sat_IDs(orbit_prop_inputs['sat_params']['data_storage_params'],sat_id_order)
            orbit_prop_inputs['sat_params']['sats_state_sorted'] = io_tools.sort_input_params_by_sat_IDs(orbit_prop_inputs['sat_params']['initial_state'],sat_id_order)
        else:
            raise NotImplementedError

        #  check that it's the right version
        if not gp_general_params['version'] == "0.6":
            raise NotImplementedError

        #  check that it's the right version
        if gp_instance_params['version'] == "0.6":

            gp_instance_params['planning_params']['planning_start_dt'] = tt.iso_string_to_dt ( gp_instance_params['planning_params']['planning_start'])
            gp_instance_params['planning_params']['planning_fixed_end_dt'] = tt.iso_string_to_dt ( gp_instance_params['planning_params']['planning_fixed_end'])
            gp_instance_params['planning_params']['planning_end_obs_dt'] = tt.iso_string_to_dt ( gp_instance_params['planning_params']['planning_end_obs'])
            gp_instance_params['planning_params']['planning_end_xlnk_dt'] = tt.iso_string_to_dt ( gp_instance_params['planning_params']['planning_end_xlnk'])
            gp_instance_params['planning_params']['planning_end_dlnk_dt'] = tt.iso_string_to_dt ( gp_instance_params['planning_params']['planning_end_dlnk'])

            dlnk_end = gp_instance_params['planning_params']['planning_end_dlnk_dt']
            obs_end = gp_instance_params['planning_params']['planning_end_obs_dt']
            xlnk_end = gp_instance_params['planning_params']['planning_end_xlnk_dt']
            if not (dlnk_end >= obs_end and dlnk_end >= xlnk_end):
                raise RuntimeWarning('Planning window end for dlnk (%s) should be set equal or later than end for observations, crosslinks (%s)'%(dlnk_end,obs_xlnk_end))

            if gp_instance_params["sats_state"] is None:
                gp_instance_params["sats_state_sorted"] = orbit_prop_inputs['sat_params']['sats_state_sorted']
            else:
                gp_instance_params["sats_state_sorted"] = io_tools.sort_input_params_by_sat_IDs(gp_instance_params["sats_state"],sat_id_order)
        else:
            raise NotImplementedError

        #  check that it's the right version
        if not data_rates_inputs['version'] == "0.3":
            raise NotImplementedError

        gp_params['orbit_prop_params'] = orbit_prop_params
        gp_params['orbit_link_params'] = orbit_link_params
        gp_params['gp_general_params'] = gp_general_params
        gp_params['gp_instance_params'] = gp_instance_params
        gp_params['data_rates_params'] = data_rates_params
        gp_params['gp_other_params'] = gp_other_params
        gp_runner = GlobalPlannerRunner (gp_params)

        # get data related to existing routes, if they were provided
        existing_route_data = data.get('existing_route_data',{})
        scheduled_routes,all_updated_routes,viz_outputs,latest_dr_uid = gp_runner.run (existing_route_data,verbose)

        output = {}
        output['version'] = OUTPUT_JSON_VER
        output['scenario_params'] = data['orbit_prop_inputs']['scenario_params']
        output['viz_data'] = viz_outputs
        output['scheduled_routes'] = scheduled_routes
        output['all_updated_routes'] = all_updated_routes
        output['latest_dr_uid'] = latest_dr_uid
        output['update_time'] = datetime.utcnow().isoformat()

        return output


def remote_multiproc_run(self,mp_queue,verbose=False):
    """ Runs the PipelineRunner as a remote multiprocessing call"""

    pr = PipelineRunner()
    print_verbose('runner_gp: remote_multiproc_run()',verbose)
    # sys.stdout = open(str(os.getpid()) + ".out", "w")
    # sys.stderr = open(str(os.getpid()) + "_error.out", "w")
    # print('stdout initialized')
    inputs = mp_queue.get()
    output = pr.run(inputs,verbose)
    mp_queue.put(output)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description='global planner')
    ap.add_argument('--prop_inputs_file',
                    type=str,
                    default='orbit_prop_inputs.json',
                    help='specify orbit propagation inputs file')

    ap.add_argument('--data_rates_file',
                    type=str,
                    default='data_rates_output.json',
                    help='specify data rate outputs file from orbit link repo')

    ap.add_argument('--link_inputs_file',
                    type=str,
                    default='orbit_link_inputs_ex.json',
                    help='specify orbit link inputs file from orbit link repo')

    ap.add_argument('--gp_general_inputs_file',
                    type=str,
                    default='crux/config/examples/gp_general_params_inputs_ex.json',
                    help='specify global planner general params file')

    ap.add_argument('--gp_inst_inputs_file',
                    type=str,
                    default=None,
                    help='specify global planner instance params file')

    ap.add_argument('--rs_s1_pickle',
                    type=str,
                    default=None,
                    help='specify post route selection step1 pickle to load')

    ap.add_argument('--rs_s2_pickle',
                    type=str,
                    default=None,
                    help='specify post route selection step2 pickle to load')

    ap.add_argument('--as_pickle',
                    type=str,
                    default=None,
                    help='specify post activity scheduling pickle to load')

    #  added this to deal with trailing white space at the command line.  freaking argparse can't understand the fact that two trailing spaces is not another argument....
    # ap.add_argument('args', nargs=argparse.REMAINDER)

    args = ap.parse_args()

    # print(args)

    pr = PipelineRunner()

    if sys.platform == 'win32':
        # todo: this is probably not the right way to support windows users...
        args.prop_inputs_file = args.prop_inputs_file.replace('/','\\')
        args.data_rates_file = args.data_rates_file.replace('/','\\')
        args.link_inputs_file = args.link_inputs_file.replace('/','\\')
        args.gp_general_inputs_file = args.gp_general_inputs_file.replace('/','\\')
        args.gp_inst_inputs_file = args.gp_inst_inputs_file.replace('/','\\')
        args.rs_s1_pickle = args.rs_s1_pickle.replace('/','\\')
        args.rs_s2_pickl = args.rs_s2_pickle.replace('/','\\')
        args.as_pickle = args.as_pickle.replace('/','\\')

    # with open(os.path.join(REPO_BASE,'crux/config/examples/orbit_prop_inputs_ex.json'),'r') as f:
    with open(os.path.join(REPO_BASE, args.prop_inputs_file),'r') as f:
        orbit_prop_inputs = json.load(f)

    with open(os.path.join(REPO_BASE,args.data_rates_file),'r') as f:
        data_rates_inputs = json.load(f)

    with open(os.path.join(REPO_BASE, args.link_inputs_file),'r') as f:
        orbit_link_inputs = json.load(f)

    with open(os.path.join(REPO_BASE,args.gp_general_inputs_file),'r') as f:
        gp_general_params_inputs = json.load(f)

    with open(os.path.join(REPO_BASE,args.gp_inst_inputs_file),'r') as f:
        gp_instance_params_inputs = json.load(f)

    data = {
        # "orbit_prop_data": orbit_prop_data,
        "orbit_prop_inputs": orbit_prop_inputs,
        "orbit_link_inputs": orbit_link_inputs,
        "gp_general_params_inputs": gp_general_params_inputs,
        "gp_instance_params_inputs": gp_instance_params_inputs,
        # "viz_params": viz_params,
        "data_rates_inputs": data_rates_inputs,
        "rs_s1_pickle": args.rs_s1_pickle if args.rs_s1_pickle != '' else None,
        "rs_s2_pickle": args.rs_s2_pickle if args.rs_s2_pickle != '' else None,
        "as_pickle": args.as_pickle if args.as_pickle != '' else None,
        "file_params":  {'new_pickle_file_name_pre': args.prop_inputs_file.split('/')[-1].split ('.')[0]}
    }

    a = time.time()
    output = pr.run(data,verbose=True)
    b = time.time()

    with open('gp_outputs.json','w') as f:
        # todo: can't currently jsonify routing objects
        del output['scheduled_routes']
        del output['all_updated_routes']
        json.dump(output ,f)

    print('run time: %f'%(b-a))