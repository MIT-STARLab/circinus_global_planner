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
import pickle
import numpy as np
from copy import deepcopy
from collections import OrderedDict
import multiprocessing as mp


#  local repo includes. todo:  make this less hackey
sys.path.append ('..')
from circinus_tools  import time_tools as tt
from circinus_tools  import io_tools
from gp_tools.io_processing import GPProcessorIO
from gp_tools.gp_plotting import GPPlotting
import gp_tools.gp_route_selection_v1 as gprsv1
import gp_tools.gp_route_selection_v2 as gprsv2
from gp_tools.gp_activity_scheduling import GPActivityScheduling
from gp_tools.gp_metrics import GPMetrics
from gp_tools.network_sim.gp_network_sim import GPNetSim
from gp_tools.routing_objects import DataMultiRoute
import gp_tools.gp_rs_execution as gp_rs_exec

# TODO: remove this line if not needed
# from gp_tools.custom_activity_window import ObsWindow

REPO_BASE = os.path.abspath(os.pardir)  # os.pardir aka '..'

OUTPUT_JSON_VER = '0.1'


class GlobalPlannerRunner:
    """easy interface for running the global planner scheduling algorithm"""

    def __init__(self, gp_params):
        """initializes based on parameters

        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """

        self.params = gp_params
        self.scenario_params = self.params['gp_orbit_prop_params']['scenario_params']
        self.sat_params = self.params['gp_orbit_prop_params']['sat_params']
        self.obs_params = self.params['gp_orbit_prop_params']['obs_params']
        self.sat_orbit_params = self.params['gp_orbit_prop_params']['sat_orbit_params']
        self.obs_params = self.params['gp_orbit_prop_params']['obs_params']
        self.pickle_params = self.params['gp_general_params']['pickle_params']
        self.other_params = self.params['gp_other_params']
        self.plot_params = self.params['gp_general_params']['plot_params']
        self.general_other_params = self.params['gp_general_params']['other_params']
        self.rs_general_params = gp_params['gp_general_params']['route_selection_general_params']
        self.rs_v2_params = gp_params['gp_general_params']['route_selection_params_v2']
        self.as_params = self.params['gp_general_params']['activity_scheduling_params']
        self.as_inst_params = self.params['gp_instance_params']['activity_scheduling_params']
        self.io_proc =GPProcessorIO(self.params)
        self.gp_plot =GPPlotting( self.params)

        if self.as_inst_params['start_utc'] < self.scenario_params['start_utc']:
            raise RuntimeError("Activity scheduling instance start time (%s) is less than scenario start time (%s)"%(self.as_inst_params['start_utc'],self.scenario_params['start_utc']))
        if self.as_inst_params['end_utc'] > self.scenario_params['end_utc']:
            raise RuntimeError("Activity scheduling instance end time (%s) is greater than scenario end time (%s)"%(self.as_inst_params['end_utc'],self.scenario_params['end_utc']))

    def pickle_rtsel_s1_stuff(self,routes_by_obs,all_stats,route_times_s,obs_indx,obs_winds,dlnk_winds_flat,ecl_winds,window_uid):

        pickle_stuff =  {}
        pickle_stuff['routes_by_obs'] = routes_by_obs
        pickle_stuff['all_stats'] = all_stats
        pickle_stuff['route_times_s'] = route_times_s
        pickle_stuff['params'] =  self.params
        pickle_stuff['obs_indx'] = obs_indx
        pickle_stuff['obs_winds'] = obs_winds
        pickle_stuff['dlnk_winds_flat'] = dlnk_winds_flat
        pickle_stuff['ecl_winds'] = ecl_winds
        pickle_stuff['window_uid'] = window_uid
        pickle_name ='pickles/rs_s1_%s_oi%d_%s' %( self.other_params['new_pickle_file_name_pre'],obs_indx,datetime.utcnow().isoformat().replace (':','_'))
        with open('%s.pkl' % ( pickle_name),'wb') as f:
            pickle.dump(  pickle_stuff,f)

    def unpickle_rtsel_s1_stuff( self):
        rs_s1_pickle = self.other_params['rs_s1_pickle_input']
        # rs_s1_pickle = rs_s1_pickle if rs_s1_pickle else self.pickle_params['route_selection_step1_pickle']

        p = pickle.load (open ( rs_s1_pickle,'rb'))
        #  TODO:  uncommon this  if want to reload parameters from file
        # self.params = p['params']

        return p['routes_by_obs'],p['all_stats'],p['route_times_s'],p['obs_indx'],p['obs_winds'],p['dlnk_winds_flat'],p['ecl_winds'],p['window_uid']

    def pickle_rtsel_s2_stuff(self,xlnk_winds_flat,sel_routes_by_obs,ecl_winds,obs_winds,dlnk_winds_flat,window_uid):

        pickle_stuff =  {}
        pickle_stuff['xlnk_winds_flat'] = xlnk_winds_flat
        pickle_stuff['sel_routes_by_obs'] = sel_routes_by_obs
        pickle_stuff['ecl_winds'] = ecl_winds
        pickle_stuff['obs_winds'] = obs_winds
        pickle_stuff['dlnk_winds_flat'] = dlnk_winds_flat
        pickle_stuff['window_uid'] = window_uid
        pickle_stuff['params'] =  self.params
        pickle_name ='pickles/rs_s2_%s_%s' %( self.other_params['new_pickle_file_name_pre'],datetime.utcnow().isoformat().replace (':','_'))
        with open('%s.pkl' % ( pickle_name),'wb') as f:
            pickle.dump(  pickle_stuff,f)

    def unpickle_rtsel_s2_stuff( self):
        rs_s2_pickle = self.other_params['rs_s2_pickle_input']
        # rs_s2_pickle = rs_s2_pickle if rs_s2_pickle else self.pickle_params['route_selection_step2_pickle']

        p = pickle.load (open ( rs_s2_pickle,'rb'))
        #  TODO:  uncommon this  if want to reload parameters from file
        # self.params = p['params']

        return p['xlnk_winds_flat'],p['sel_routes_by_obs'],p['ecl_winds'],p['obs_winds'],p['dlnk_winds_flat'],p['window_uid']


    def pickle_actsc_stuff(self,routes_by_obs,ecl_winds,scheduled_routes,energy_usage,data_usage,window_uid):

        pickle_stuff =  {}
        pickle_stuff['routes_by_obs'] = routes_by_obs
        pickle_stuff['ecl_winds'] = ecl_winds
        pickle_stuff['scheduled_routes'] = scheduled_routes
        pickle_stuff['energy_usage'] = energy_usage
        pickle_stuff['data_usage'] = data_usage
        pickle_stuff['window_uid'] = window_uid
        pickle_name ='pickles/as_%s' %(datetime.utcnow().isoformat().replace (':','_'))
        with open('%s.pkl' % ( pickle_name),'wb') as f:
            pickle.dump(  pickle_stuff,f)


    def unpickle_actsc_stuff( self):
        as_pickle = self.other_params['as_pickle_input']

        p = pickle.load (open ( as_pickle,'rb'))

        return p['routes_by_obs'],p['ecl_winds'],p['scheduled_routes'],p['energy_usage'],p['data_usage'],p['window_uid']



    def run_nominal_activity_scheduling( self, routes_by_obs,ecl_winds):
        gp_as = GPActivityScheduling ( self.params)

        # flatten the list of all routes, which currently has nested lists for each observation
        routes_flat = [rt for rts in routes_by_obs.values() for rt in rts]

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        print('make activity scheduling model')
        gp_as.make_model (routes_flat, ecl_winds,verbose = True)
        stats =gp_as.get_stats (verbose = True)
        print('solve activity scheduling')
        t_a = time.time()
        gp_as.solve ()
        t_b = time.time()
        # gp_as.print_sol ()
        print('extract_routes')
        #  make a copy of the windows in the extracted routes so we don't mess with the original objects ( just to be extra careful)
        routes = gp_as.extract_utilized_routes ( copy_routes = True, verbose  = False)
        energy_usage,data_usage = gp_as.extract_resource_usage(  decimation_factor =1)

        time_elapsed = t_b-t_a

        return  routes, energy_usage, data_usage

    def calc_activity_scheduling_results ( self,obs_winds,dlnk_winds_flat,rs_routes_by_obs,sched_routes, energy_usage):
        gp_met = GPMetrics(self.params)

        def in_sched_window(wind):
            return wind.start >= self.as_inst_params['start_utc_dt'] and wind.end <= self.as_inst_params['end_utc_dt']

        total_collectible_DV_all_obs_winds = sum(obs.data_vol for winds in obs_winds for obs in winds  if in_sched_window(obs))
        total_dlnkable_DV_all_dlnk_winds = sum(dlnk.data_vol for winds in dlnk_winds_flat for dlnk in winds if in_sched_window(dlnk))
        rs_output_routes = [rt for rts in rs_routes_by_obs.values() for rt in rts]
        total_throughput_DV_rs_routes = sum(sum(rt.data_vol for rt in rts) for obs, rts in rs_routes_by_obs.items() if in_sched_window(obs))
        total_collectible_DV_rs_routes = sum(min(obs.data_vol,sum(rt.data_vol for rt in rts)) for obs, rts in rs_routes_by_obs.items() if in_sched_window(obs))

        print('------------------------------')
        print('in scheduling window:')
        print('len(rs_output_routes)')
        print(len(rs_output_routes))
        print('len(sched_routes)')
        print(len(sched_routes))
        print('total_collectible_DV_all_obs_winds')
        print(total_collectible_DV_all_obs_winds)
        print('total_dlnkable_DV_all_dlnk_winds')
        print(total_dlnkable_DV_all_dlnk_winds)
        print('total_throughput_DV_rs_routes')
        print(total_throughput_DV_rs_routes)
        print('total_collectible_DV_rs_routes')
        print(total_collectible_DV_rs_routes)
        print('weights')
        print(self.as_params['obj_weights'])
        # dv_stats = gp_met.assess_dv_all_routes (sched_routes,verbose = True)
        dv_obs_stats = gp_met.assess_dv_by_obs (rs_routes_by_obs,sched_routes,verbose = True)
        lat_stats = gp_met.assess_latency_all_routes (sched_routes,verbose = True)
        lat_obs_stats = gp_met.assess_latency_by_obs (rs_routes_by_obs,sched_routes,verbose = True)
        aoi_targ_stats = gp_met.assess_aoi_by_obs_target(rs_routes_by_obs,sched_routes,verbose = True)

        gp_netsim = GPNetSim ( self.params, self.io_proc)
        gp_netsim.sim_tlm_cmd_routing(sched_routes, verbose =  False)
        #  this is indexed by sat index
        sats_cmd_update_hist = gp_netsim.get_all_sats_cmd_update_hist()
        aoi_sat_cmd_stats = gp_met.assess_aoi_sat_cmd(sats_cmd_update_hist,verbose = True)
        #  this is  indexed by ground station index
        sats_tlm_update_hist = gp_netsim.get_all_sats_tlm_update_hist()
        aoi_sat_tlm_stats = gp_met.assess_aoi_sat_tlm(sats_tlm_update_hist,verbose = True)

        resource_margin_stats = gp_met.assess_resource_margin(energy_usage,verbose = True)


        plot_outputs = {}
        plot_outputs['rs_targIDs_found'] = aoi_targ_stats['rs_targIDs_found']
        plot_outputs['obs_aoi_curves_by_targID'] = aoi_targ_stats['aoi_curves_by_targID_sched']
        plot_outputs['cmd_aoi_curves_by_sat_indx'] = aoi_sat_cmd_stats['aoi_curves_by_sat_indx']
        plot_outputs['tlm_aoi_curves_by_sat_indx'] = aoi_sat_tlm_stats['aoi_curves_by_sat_indx']
        return plot_outputs

    def run_route_selection(self,gp_ps,obs,dlnk_winds_flat,xlnk_winds,obj_weights,dr_uid):
        """[summary]

        [description]
        :param gp_ps: [description]
        :type gp_ps: [type]
        :param obs: [description]
        :type obs: [type]
        :param dlnk_winds_flat: [description]
        :type dlnk_winds_flat: [type]
        :param xlnk_winds: [description]
        :type xlnk_winds: [type]
        :param obj_weights: [description]
        :type obj_weights: [type]
        :param dr_uid: globally unique data route ID across all data routes for all observations
        :type dr_uid: int
        :returns: [description]
        :rtype: {[type]}
        """



        return routes,obs,stats,time_elapsed,dr_uid

    def run_nominal_route_selection_v1( self,obs_winds,dlnk_winds_flat,xlnk_winds):
        gp_ps = gprsv1.GPDataRouteSelection ( self.params)

        print ('nominal route selection v1')

        obj_weights = {
            "total_dv": 0.45,
            "num_paths_sel": 0.1,
            "latency_sf": 0.45,
        }

        obs_indx =0
        #  list of all routes, indexed by obs_indx ( important for pickling)
        all_routes =[]
        all_routes_obs =[]
        all_stats =[]
        route_times_s =[]

        # this should be unique across all data routes
        dr_uid = 0

        for sat_indx in range( self.sat_params['num_sats']):
            for  index, obs in  enumerate ( obs_winds[sat_indx]):

                print ("sat_indx")
                print (sat_indx)
                print ("obs")
                print ( index)

                # routes,obs,stats,time_elapsed,dr_uid = self.run_route_selection(gp_ps,obs,dlnk_winds_flat,xlnk_winds,obj_weights,dr_uid)

                # run the route selection algorithm
                gp_ps.make_model (obs,dlnk_winds_flat,xlnk_winds, obj_weights, verbose = True)
                stats =gp_ps.get_stats (verbose = True)
                t_a = time.time()
                gp_ps.solve ()
                t_b = time.time()
                gp_ps.print_sol ()
                adjust_opt = self.rs_general_params['adjust_window_timing_post_selection']
                routes,dr_uid = gp_ps. extract_routes ( dr_uid,copy_windows= True,adjust_window_timing=adjust_opt,verbose  = True)

                time_elapsed = t_b-t_a
                # end algorithm

                print ('obs.data_vol, total dlnk dv, ratio')
                print (obs.data_vol,sum(dr.data_vol for dr in routes),sum(dr.data_vol for dr in routes)/obs.data_vol)
                print ('min latency, ave latency, max latency')
                latencies = [(rt.route[-1].end - rt.route[0].end).total_seconds()/60 for rt in  routes]
                if len(latencies) > 0:
                    print (np.min( latencies),np.mean( latencies),np.max( latencies))
                else:
                    print('no routes found')
                print ('time_elapsed')
                print (time_elapsed)

                all_routes.append ( routes)
                all_routes_obs.append ( obs)
                all_stats.append ( stats)
                route_times_s.append ( time_elapsed)

                obs_indx +=1

                if obs_indx >= 1:
                    break

            if obs_indx >= 1:
                break

        return all_routes,all_routes_obs,all_stats,route_times_s,obs_indx

    def run_nominal_route_selection_v2_step1( self,obs_winds,dlnk_winds_flat,xlnk_winds,verbose=False):

        if verbose:
            print ('nominal route selection v2 step 1')

        obs_indx =0
        # dict of all routes, with obs as key
        routes_by_obs ={}
        indices_by_obs ={}
        all_stats =[]
        route_times_s =[]

        def filter_obs(obs_winds):
            obs_winds_filt = [[] for sat_indx in range( self.sat_params['num_sats'])]

            for sat_indx in range( self.sat_params['num_sats']):
                for sat_obs_index,obs_wind in enumerate(obs_winds[sat_indx]):
                    if obs_wind.start >= self.as_inst_params['start_utc_dt'] and obs_wind.end <= self.as_inst_params['end_utc_dt']:
                        obs_winds_filt[sat_indx].append(obs_wind)
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

                    # print ('len(routes)')
                    # print (len(routes))
                    # print ('obs.data_vol, total dlnk dv, ratio')
                    # print (obs.data_vol,sum(dr.data_vol for dr in routes),sum(dr.data_vol for dr in routes)/obs.data_vol)
                    # print ('min latency, ave latency, max latency')
                    # latencies = [(rt.route[-1].end - rt.route[0].end).total_seconds()/60 for rt in  routes]
                    # if len(latencies) > 0:
                    #     print (np.min( latencies),np.mean( latencies),np.max( latencies))
                    # else:
                    #     print('no routes found')

        t_b = time.time()
        time_elapsed = t_b-t_a

        # this should be unique across all data routes
        dr_uid = 0

        # explicitly validate routes
        for routes in routes_by_obs.values():
            for dr in routes:
                #  set the data route ID now, because they were not set uniquely if we were running in parallel
                dr.ID = dr_uid
                dr_uid += 1

                # try:
                dr.validate(dv_epsilon = self.as_params['dv_epsilon_Mb'])
                # except:
                #     raise Exception("Couldn't handle dr %s for obs %s (indices %s)"%(dr,dr.get_obs(),indices_by_obs[dr.get_obs()]))

        if verbose:
            print('total_rs_step1_runtime')
            print(time_elapsed)
            print('runtime per rs iteration')
            print('ave: %fs'%(np.mean(route_times_s)))
            print('std: %fs'%(np.std(route_times_s)))

        return routes_by_obs,all_stats,route_times_s,obs_indx, dr_uid

    def  setup_test( self,obs_winds,dlnk_winds_flat,xlnk_winds):
        obs_winds_sel = [[] for sat_indx in range (self.sat_params['num_sats'])]
        # obs_winds_sel[13].append ( obs_winds[13][0])
        # obs_winds_sel[19].append ( obs_winds[19][0])
        # obs_winds_sel[19].append ( obs_winds[19][1])
        # obs_winds_sel[20].append ( obs_winds[20][6])
        # obs_winds_sel[26].append ( obs_winds[26][0])

        if False:
            self.gp_plot.plot_winds(
                range (self.sat_params['num_sats']),
                # [13,19,20,26],
                obs_winds_sel,
                [],
                dlnk_winds_flat,
                [],
                [],
                {},
                self.scenario_params['start_utc_dt'],
                self.scenario_params['end_utc_dt'],
                plot_title = 'all obs, dlnks',
                plot_size_inches = (18,12),
                plot_include_labels = True,
                show=  False,
                fig_name='plots/temp1.pdf'
            )

        return obs_winds_sel


    def run_test_route_selection( self,obs_winds,dlnk_winds_flat,xlnk_winds):
        gp_rs = GPDataRouteSelection ( self.params)

        total_dv_weights = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        num_paths_sel_weights = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
        latency_sf_weights = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]
        weights_tups = zip(total_dv_weights,num_paths_sel_weights,latency_sf_weights)

        print ('test route selection')

        # obs_winds_sel =  self.setup_test(obs_winds,dlnk_winds_flat,xlnk_winds)
        obs_winds_sel = [[] for sat_indx in range (self.sat_params['num_sats'])]
        # obs_winds_sel[13].append ( obs_winds[13][0])
        # obs_winds_sel[19].append ( obs_winds[19][0])
        # obs_winds_sel[19].append ( obs_winds[19][1])
        obs_winds_sel[20].append ( obs_winds[20][6])
        # obs_winds_sel[26].append ( obs_winds[26][0])

        obs_indx =0
        #  list of all routes, indexed by obs_indx ( important for pickling)
        all_routes =[]
        all_routes_obs =[]
        all_stats =[]
        route_times_s =[]

        #  flatten into a simple list of observations
        obs_winds_sel_flat = flat_list = [item for sublist in obs_winds_sel for item in sublist]


        for  index, obs in  enumerate ( obs_winds_sel_flat):
            for  weights_index, weights in  enumerate(weights_tups):
                print ("obs")
                print ( index)
                print (  'weights_index')
                print (  weights_index)

                obj_weights = {
                    "total_dv":   weights[0],
                    "num_paths_sel":   weights[1],
                    "latency_sf":   weights[2],
                }

                routes,obs,stats,time_elapsed = self.run_route_selection(gp_rs,obs,dlnk_winds_flat,xlnk_winds,obj_weights)

                print ('obs.data_vol, total dlnk dv, ratio')
                print (obs.data_vol,sum(dr.data_vol for dr in routes),sum(dr.data_vol for dr in routes)/obs.data_vol)
                print ('min latency, ave latency, max latency')
                latencies = [(rt.route[-1].end - rt.route[0].end).total_seconds()/60 for rt in  routes]
                print (np.min( latencies),np.mean( latencies),np.max( latencies))
                print ('time_elapsed')
                print (time_elapsed)

                all_routes.append ( routes)
                all_routes_obs.append ( obs)
                all_stats.append ( stats)
                route_times_s.append ( time_elapsed)

            obs_indx +=1

        return all_routes,all_routes_obs,all_stats,route_times_s,obs_indx,weights_tups

    def pareto_plot(all_routes):
        total_dv_weights = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        num_paths_sel_weights = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
        latency_sf_weights = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]
        weights_tups = zip(total_dv_weights,num_paths_sel_weights,latency_sf_weights)
        self.gp_plot.plot_route_latdv_pareto(all_routes,weights_tups,'plots/obs_winds_20_6_pareto.pdf')

    def plot_route_selection_results ( self,routes_by_obs,dlnk_winds_flat,xlnk_winds_flat,num_obs_to_plot):

        first_dlnk_indx = 10
        num_dlnks_to_plot = 20
        for obs_indx, (obs,routes) in enumerate (routes_by_obs.items()):
            if obs_indx >= num_obs_to_plot:
                break

            print('Plotting obs index %d'%(obs_indx))

            rts_by_dlnk = OrderedDict()

            for dr in routes:
                dlnk = dr.get_dlnk()
                rts_by_dlnk.setdefault(dlnk,[]).append(dr)

            for dlnk_indx, dlnk in enumerate (rts_by_dlnk.keys()):
                if dlnk_indx >= num_dlnks_to_plot:
                    break
                if dlnk_indx < first_dlnk_indx:
                    continue

                # if dlnk_indx != 33:
                #     continue

                rts = rts_by_dlnk[dlnk]

                # TODO:  this stuff needs to be  changed
                # (sel_obs_winds_flat, sel_dlnk_winds_flat, sel_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind) = self.io_proc.extract_flat_windows (routes)
                (sel_obs_winds_flat, sel_dlnk_winds_flat, sel_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind) = self.io_proc.extract_flat_windows (rts)

                obs = rts[0].get_obs()

                sats_to_include =  range (self.sat_params['num_sats'])
                # sats_to_include = [1,2,3]

                plot_len = max(self.rs_general_params['wind_filter_duration_s'],self.rs_general_params['wind_filter_duration_obs_sat_s'])

                #  plot the selected down links and cross-links
                self.gp_plot.plot_winds(
                    sats_to_include,
                    sel_obs_winds_flat,
                    sel_obs_winds_flat,
                    dlnk_winds_flat,
                    sel_dlnk_winds_flat,
                    # [],
                    xlnk_winds_flat,
                    sel_xlnk_winds_flat,
                    # [],
                    route_indcs_by_wind,
                    self.scenario_params['start_utc_dt'],
                    # obs.start,
                    obs.start + timedelta( seconds= plot_len),
                    # self.scenario_params['start_utc_dt'],
                    # self.scenario_params['start_utc_dt'] + timedelta( seconds= self.rs_general_params['wind_filter_duration_s']),
                    # self.scenario_params['end_utc_dt']-timedelta(minutes=200),
                    plot_title = 'Route Plot',
                    plot_size_inches = (18,12),
                    plot_include_dlnk_labels = self.rs_general_params['plot_include_dlnk_labels'],
                    plot_include_xlnk_labels = self.rs_general_params['plot_include_xlnk_labels'],
                    show= False,
                    fig_name='plots/test_rs_o{0}_d{1}_6sat.pdf'.format (obs_indx,dlnk_indx)
                )


    def  plot_activity_scheduling_results ( self,all_possible_winds,routes_by_obs,routes,energy_usage,data_usage,ecl_winds,metrics_plot_inputs):

        # do a bunch of stuff to extract the windows from all of the routes as indexed by observation
        # note that this stuff is not thewindows from the scheduled routes, but rather the windows from all the route selected in route selection
        #  start
        sel_obs_winds_flat = [set() for  sat_indx  in range  (self.sat_params['num_sats'])]
        sel_dlnk_winds_flat = [set() for sat_indx  in range (self.sat_params['num_sats'])]
        sel_xlnk_winds_flat = [set() for sat_indx  in range (self.sat_params['num_sats'])]

        for rts_indx, (obs,rts) in enumerate (routes_by_obs.items()):
            obs_winds_rt, dlnk_winds_rt, \
            xlnk_winds_rt, _, _ = self.io_proc.extract_flat_windows (rts)

            for sat_indx in range  (self.sat_params['num_sats']):
                [sel_obs_winds_flat[sat_indx] .add( wind)  for wind in obs_winds_rt[sat_indx]]
                [sel_dlnk_winds_flat[sat_indx] .add(wind ) for wind in dlnk_winds_rt[sat_indx]]
                [sel_xlnk_winds_flat[sat_indx] .add(wind ) for wind in xlnk_winds_rt[sat_indx]]

        for sat_indx in range  (self.sat_params['num_sats']):
            sel_obs_winds_flat[sat_indx] = list(sel_obs_winds_flat[sat_indx])
            sel_dlnk_winds_flat[sat_indx] = list(sel_dlnk_winds_flat[sat_indx])
            sel_xlnk_winds_flat[sat_indx] = list(sel_xlnk_winds_flat[sat_indx])
            sel_obs_winds_flat[sat_indx].sort(key=lambda x: x.start)
            sel_dlnk_winds_flat[sat_indx].sort(key=lambda x: x.start)
            sel_xlnk_winds_flat[sat_indx].sort(key=lambda x: x.start)
        # end

        sched_obs_winds_flat, sched_dlnk_winds_flat, \
        sched_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind = self.io_proc.extract_flat_windows (routes,copy_windows= False)

        #
        sats_to_include =  [sat_id for sat_id in self.sat_params['sat_id_order']]
        # sats_to_include =  [sat_id for sat_id in range(20,30)]
        # sats_to_include = [12,13,14,15,16]

        all_obs_winds,all_dlnk_winds_flat,all_xlnk_winds_flat = all_possible_winds

        # plot all winds
        self.gp_plot.plot_winds(
            sats_to_include,
            all_obs_winds,
            all_obs_winds,
            all_dlnk_winds_flat,
            all_dlnk_winds_flat,
            all_xlnk_winds_flat,
            None,
            None,
            self.as_inst_params['start_utc_dt'],
            # self.as_inst_params['start_utc_dt'] + timedelta( seconds= self.rs_general_params['wind_filter_duration_s']),
            self.as_inst_params['end_utc_dt'],
            base_time = self.scenario_params['start_utc_dt'],
            plot_title = 'All Possible Activities',
            plot_size_inches = (18,12),
            plot_include_dlnk_labels = self.as_params['plot_include_dlnk_labels'],
            plot_include_xlnk_labels = self.as_params['plot_include_xlnk_labels'],
            show=  False,
            fig_name='plots/test_all_windows.pdf'
        )


        # plot the selected down links and cross-links this
        self.gp_plot.plot_winds(
            sats_to_include,
            sel_obs_winds_flat,
            sched_obs_winds_flat,
            sel_dlnk_winds_flat,
            sched_dlnk_winds_flat,
            sel_xlnk_winds_flat,
            sched_xlnk_winds_flat,
            route_indcs_by_wind,
            self.as_inst_params['start_utc_dt'],
            # self.as_inst_params['start_utc_dt'] + timedelta( seconds= self.rs_general_params['wind_filter_duration_s']),
            self.as_inst_params['end_utc_dt'],
            base_time = self.scenario_params['start_utc_dt'],
            plot_title = 'Scheduled Activities',
            plot_size_inches = (18,12),
            plot_include_dlnk_labels = self.as_params['plot_include_dlnk_labels'],
            plot_include_xlnk_labels = True, #self.as_params['plot_include_xlnk_labels'],
            show=  False,
            fig_name='plots/test_activity_times.pdf'
        )

        # self.gp_plot.plot_data_circles(
        #     sats_to_include,
        #     sched_obs_winds_flat,
        #     sched_obs_winds_flat,
        #     sched_dlnk_winds_flat,
        #     sched_dlnk_winds_flat,
        #     sched_xlnk_winds_flat,
        #     sched_xlnk_winds_flat,
        #     route_indcs_by_wind,
        #     self.as_inst_params['start_utc_dt'],
        #     # self.as_inst_params['start_utc_dt'] + timedelta( seconds= self.rs_general_params['wind_filter_duration_s']),
        #     self.as_inst_params['end_utc_dt'],
              # base_time = self.scenario_params['start_utc_dt'],
        #     plot_title = 'Activity Data Volumes',
        #     plot_size_inches = (18,12),
        #     plot_include_labels = self.as_params['plot_include_labels'],
        #     show=  False,
        #     fig_name='plots/test_data_volume.pdf'
        # )

        self.gp_plot.plot_energy_usage(
            sats_to_include,
            energy_usage,
            ecl_winds,
            self.as_inst_params['start_utc_dt'],
            self.as_inst_params['end_utc_dt'],
            base_time = self.scenario_params['start_utc_dt'],
            plot_title = 'Energy Storage Utilization',
            plot_size_inches = (18,12),
            show=  False,
            fig_name='plots/test_energy.pdf'
        )

        # note that little blips upward that appear in data storage plot for "just passing through" crosslink pairs (looks like _____^_____ . Isn't that some great asciiart?) can be caused by a slight overlap in timepoints at beg/end of two back to back activities. This should not cause any problems
        self.gp_plot.plot_data_usage(
            sats_to_include,
            data_usage,
            ecl_winds,
            self.as_inst_params['start_utc_dt'],
            self.as_inst_params['end_utc_dt'],
            base_time = self.scenario_params['start_utc_dt'],
            plot_title = 'Data Storage Utilization',
            plot_size_inches = (18,12),
            show=  False,
            fig_name='plots/test_data.pdf'
        )

        found_targIDs = metrics_plot_inputs['rs_targIDs_found']
        # targs_to_include = [targ['id'] for targ in self.obs_params['targets']]
        # targs_to_include = [0,3,4,7,8]
        # targs_to_include = range(15)
        # targs_to_include = found_targIDs[0:10]
        targs_to_include = found_targIDs

        self.gp_plot.plot_obs_aoi(
            targs_to_include,
            metrics_plot_inputs['obs_aoi_curves_by_targID'],
            self.as_inst_params['start_utc_dt'],
            self.as_inst_params['end_utc_dt'],
            base_time = self.scenario_params['start_utc_dt'],
            plot_title = 'Observation Target AoI',
            plot_size_inches = (18,12),
            show=False,
            fig_name='plots/test_obs_aoi_plot.pdf'
        )

        # sats_to_include =  [sat_p['sat_id'] for sat_p in self.sat_orbit_params]
        # sats_to_include =  range(10)
        aoi_option = 'cmd'
        self.gp_plot.plot_sat_tlm_cmd_aoi(
            sats_to_include,
            metrics_plot_inputs['cmd_aoi_curves_by_sat_indx'],
            aoi_option,
            self.as_inst_params['start_utc_dt'],
            self.as_inst_params['end_utc_dt'],
            base_time = self.scenario_params['start_utc_dt'],
            plot_title = 'Satellite Command Uplink AoI',
            plot_size_inches = (18,12),
            show=False,
            fig_name='plots/test_cmd_aoi_plot.pdf'
        )

        aoi_option = 'tlm'
        self.gp_plot.plot_sat_tlm_cmd_aoi(
            sats_to_include,
            metrics_plot_inputs['tlm_aoi_curves_by_sat_indx'],
            aoi_option,
            self.as_inst_params['start_utc_dt'],
            self.as_inst_params['end_utc_dt'],
            base_time = self.scenario_params['start_utc_dt'],
            plot_title = 'Satellite Telemetry Downlink AoI',
            plot_size_inches = (18,12),
            show=False,
            fig_name='plots/test_tlm_aoi_plot.pdf'
        )


    def validate_unique_windows( self,obs_winds,dlnk_winds_flat,xlnk_winds,ecl_winds):
        all_wind_ids = set()

        for indx in range(len(obs_winds)):
            for wind in obs_winds[indx]:
                if wind.window_ID in all_wind_ids:
                    raise Exception('Found a duplicate unique window ID where it should not have been possible')
                all_wind_ids.add(wind.window_ID)

        for indx in range(len(dlnk_winds_flat)):
            for wind in dlnk_winds_flat[indx]:
                if wind.window_ID in all_wind_ids:
                    raise Exception('Found a duplicate unique window ID where it should not have been possible')
                all_wind_ids.add(wind.window_ID)

        for indx in range(len(xlnk_winds)):
            for indx_2 in range(len(xlnk_winds[indx])):
                for wind in xlnk_winds[indx][indx_2]:
                    if wind.window_ID in all_wind_ids:
                        raise Exception('Found a duplicate unique window ID where it should not have been possible')
                    all_wind_ids.add(wind.window_ID)

        for indx in range(len(ecl_winds)):
            for wind in ecl_winds[indx]:
                if wind.window_ID in all_wind_ids:
                    raise Exception('Found a duplicate unique window ID where it should not have been possible')
                all_wind_ids.add(wind.window_ID)

    def run_nominal_route_selection_v2_step2(self,routes_by_obs):
        print('num routes')
        print(sum(len(rts) for rts in routes_by_obs.values()))


        gp_met = GPMetrics(self.params)
        # t_a = time.time()
        # need to figure out how to window this or something so that we don't have to compare to every other data route - that's horribly expensive
        print('Assess route overlap pre RS step 2')
        overlap_cnt_by_route,stats_rs1 = gp_met.assess_route_overlap( routes_by_obs,verbose=True)
        # t_b = time.time()
        # time_elapsed = t_b-t_a
        # print('time_elapsed')
        # print(time_elapsed)

        gp_rs = gprsv2.GPDataRouteSelection ( self.params)

        selected_rts_by_obs = gp_rs.run_step2(routes_by_obs,overlap_cnt_by_route)

        print('Assess route overlap post RS step 2')
        overlap_cnt_by_route,stats_rs2 = gp_met.assess_route_overlap( selected_rts_by_obs,verbose=True)


        for rts in selected_rts_by_obs.values():
            for dmr in rts:
                dmr.validate(self.as_params['dv_epsilon_Mb'])

        return selected_rts_by_obs

    def run( self):

        #################################
        #  parse inputs, if desired
        #################################

        if self.general_other_params['load_windows_from_file']:
            print('Load files')

            # parse the inputs into activity windows
            window_uid = 0
            print('Load obs')
            obs_winds, window_uid =self.io_proc.import_obs_winds(window_uid)
            print('Load dlnks')
            dlnk_winds, dlnk_winds_flat, window_uid =self.io_proc.import_dlnk_winds(window_uid)
            print('Load xlnks')
            xlnk_winds, xlnk_winds_flat, window_uid =self.io_proc.import_xlnk_winds(window_uid)
            print('Load ecl')
            ecl_winds, window_uid =self.io_proc.import_eclipse_winds(window_uid)

            # with open('temp.pkl','wb') as f:
            #     pickle.dump( {'params': self.params},f)

            # important to check this because window unique IDs are used as hashes in dictionaries in the scheduling code
            print('Validate windows')
            self.validate_unique_windows(obs_winds,dlnk_winds_flat,xlnk_winds,ecl_winds)

            # todo:  probably ought to delete the input times and rates matrices to free up space

            print('In windows loaded from file:')
            print('obs_winds')
            print(sum([len(p) for p in obs_winds]))
            print('dlnk_win')
            print(sum([len(p) for p in dlnk_winds]))
            print('xlnk_win')
            print(sum([len(xlnk_winds[i][j]) for i in  range( self.sat_params['num_sats']) for j in  range( self.sat_params['num_sats']) ]))

        #################################
        #  route selection step 1
        #################################

        pas_a = time.time()

        # If we need output from step 1
        run_step_1 = not self.other_params['rs_s2_pickle_input'] and not self.other_params['as_pickle_input']
        if run_step_1:

            #  if  we are loading from file, do that
            if self.other_params['rs_s1_pickle_input']:
                print('Unpickling route selection step one stuff')
                routes_by_obs,all_stats,route_times_s,obs_indx,obs_winds,dlnk_winds_flat,ecl_winds,window_uid = self.unpickle_rtsel_s1_stuff()

            #  otherwise run route selection step 1
            else:
                print('Run route selection step 1')
                routes_by_obs,all_stats,route_times_s, obs_indx, dr_uid  =  self.run_nominal_route_selection_v2_step1(obs_winds,dlnk_winds_flat,xlnk_winds,verbose=self.rs_v2_params['verbose_step1'])
                # routes_by_obs,all_stats,route_times_s, obs_indx, weights_tups  =  self.run_test_route_selection(obs_winds,dlnk_winds_flat,xlnk_winds)

            #  pickle before step 2 because step 2 doesn't take that long
            if self.pickle_params['pickle_route_selection_step1_results']:
                self.pickle_rtsel_s1_stuff(routes_by_obs,all_stats,route_times_s,obs_indx,obs_winds,dlnk_winds_flat,ecl_winds,window_uid)
        else:
            print('Skipping route selection step one stuff')


        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        #################################
        #  route selection step 2
        #################################

        # If we need output from step 2
        run_step_2 = not self.other_params['as_pickle_input']
        if run_step_2:
            if self.other_params['rs_s2_pickle_input']:
                print('Unpickling route selection step two stuff')
                xlnk_winds_flat,sel_routes_by_obs,ecl_winds,obs_winds,dlnk_winds_flat,window_uid = self.unpickle_rtsel_s2_stuff()
            else:
                print('np.mean(route_times_s)')
                print(np.mean(route_times_s))
                print('np.std(route_times_s)')
                print(np.std(route_times_s))
                passthru = False
                if passthru:
                    sel_routes_by_obs = {obs:[DataMultiRoute(ID=0,data_routes=[dr]) for dr in rts] for obs,rts in routes_by_obs.items()}
                else:
                    # run step 2. todo:  move this elsewhere
                    sel_routes_by_obs = self.run_nominal_route_selection_v2_step2(routes_by_obs)

            if self.pickle_params['pickle_route_selection_step2_results']:
                self.pickle_rtsel_s2_stuff(xlnk_winds_flat,sel_routes_by_obs,ecl_winds,obs_winds,dlnk_winds_flat,window_uid)
        else:
            print('Skipping route selection step two stuff')

        #################################
        # route selection output stage
        #################################

        print('route selection output stage')

        if self.rs_general_params['plot_route_selection_results']:
            self.plot_route_selection_results (sel_routes_by_obs,dlnk_winds_flat,xlnk_winds_flat,num_obs_to_plot = 5)

        #################################
        #  activity scheduling stage
        #################################

        if not self.as_params['run_activity_scheduling']:
            return None

        if self.other_params['as_pickle_input']:
            sel_routes_by_obs,ecl_winds,scheduled_routes,energy_usage,data_usage, window_uid = self.unpickle_actsc_stuff()
        else:
            found_routes = any([len(rts) >0 for rts in sel_routes_by_obs.values()])

            #  to protect against the weird case where we didn't find any routes ( shouldn't happen, unless we're at the very end of the simulation, or you're trying to break things)
            if found_routes:
                scheduled_routes,energy_usage,data_usage = self.run_nominal_activity_scheduling(sel_routes_by_obs,ecl_winds)
            else:
                scheduled_routes,energy_usage,data_usage = ([],None,None)
                print('No routes were found in route selection; not running activity selection')

        # if we are saving to file, do that
        if self.pickle_params['pickle_act_scheduling_results']:
            self.pickle_actsc_stuff(sel_routes_by_obs,ecl_winds,scheduled_routes,energy_usage,data_usage, window_uid)


        pas_b = time.time()
        total_plan_and_sched_runtime = pas_b - pas_a

        metrics_plot_inputs = self.calc_activity_scheduling_results (obs_winds,dlnk_winds_flat,sel_routes_by_obs,scheduled_routes, energy_usage)

        print('total_plan_and_sched_runtime (warning: may include (un)pickling time and RS plot output)')
        print("%.2f seconds"%(total_plan_and_sched_runtime))

        #################################
        #  Activity scheduling output stage
        #################################

        print('activity scheduling output stage')

        if self.as_params['plot_activity_scheduling_results']:
            self.plot_activity_scheduling_results((obs_winds,dlnk_winds_flat,xlnk_winds_flat),sel_routes_by_obs,scheduled_routes,energy_usage,data_usage,ecl_winds,metrics_plot_inputs)


        # if you want to see windows from RS output...
        # sel_routes_flat = [dr for rts in sel_routes_by_obs.values() for dr in rts]
        # (sel_obs_winds_flat, sel_dlnk_winds_flat, sel_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind) = self.io_proc.extract_flat_windows (sel_routes_flat)
        # outputs= self.io_proc.make_sat_history_outputs (sel_obs_winds_flat, sel_xlnk_winds_flat, sel_dlnk_winds_flat, link_info_by_wind)

        (sched_obs_winds_flat, sched_dlnk_winds_flat, sched_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind) = self.io_proc.extract_flat_windows (scheduled_routes)
        outputs= self.io_proc.make_sat_history_outputs (sched_obs_winds_flat, sched_xlnk_winds_flat, sched_dlnk_winds_flat, link_info_by_wind)


        return outputs


class PipelineRunner:

    def run(self, data):
        """

        """

        orbit_prop_inputs = deepcopy( data['orbit_prop_inputs'])
        orbit_link_inputs = data['orbit_link_inputs']
        gp_general_params_inputs = data['gp_general_params_inputs']
        gp_instance_params_inputs = data['gp_instance_params_inputs']
        data_rates_inputs = data['data_rates_inputs']
        file_params = data.get ('file_params',{})

        gp_params = {}
        gp_orbit_prop_params = orbit_prop_inputs
        gp_orbit_link_params = orbit_link_inputs
        gp_general_params = gp_general_params_inputs
        gp_instance_params = gp_instance_params_inputs
        gp_data_rates_params = data_rates_inputs
        gp_other_params = {}
        gp_other_params['new_pickle_file_name_pre']  = file_params.get ('orbit_prop_inputs_file' ,'default.json').split ('.')[0]
        gp_other_params['rs_s1_pickle_input']  = data['rs_s1_pickle']
        gp_other_params['rs_s2_pickle_input']  = data['rs_s2_pickle']
        gp_other_params['as_pickle_input']  = data['as_pickle']

        if data['rs_s1_pickle'] and data['rs_s2_pickle']:
            raise Exception('Should only specify 1 input pickle for route selection')

        if orbit_prop_inputs['version'] == "0.4":
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
            orbit_prop_inputs['sat_params']['initial_state_sorted'] = io_tools.sort_input_params_by_sat_IDs(orbit_prop_inputs['sat_params']['initial_state'],sat_id_order)
        else:
            raise NotImplementedError

        #  check that it's the right version
        if not gp_general_params['version'] == "0.1":
            raise NotImplementedError

        #  check that it's the right version
        if gp_instance_params['version'] == "0.1":
            # gp_instance_params['route_selection_params']['start_utc_dt'] = tt.iso_string_to_dt ( gp_instance_params['route_selection_params']['start_utc'])
            gp_instance_params['activity_scheduling_params']['start_utc_dt'] = tt.iso_string_to_dt ( gp_instance_params['activity_scheduling_params']['start_utc'])
            gp_instance_params['metrics_params']['start_utc_dt'] = tt.iso_string_to_dt ( gp_instance_params['metrics_params']['start_utc'])
            # gp_instance_params['route_selection_params']['end_utc_dt'] = tt.iso_string_to_dt ( gp_instance_params['route_selection_params']['end_utc'])
            gp_instance_params['activity_scheduling_params']['end_utc_dt'] = tt.iso_string_to_dt ( gp_instance_params['activity_scheduling_params']['end_utc'])
            gp_instance_params['metrics_params']['end_utc_dt'] = tt.iso_string_to_dt ( gp_instance_params['metrics_params']['end_utc'])
        else:
            raise NotImplementedError

        #  check that it's the right version
        if not data_rates_inputs['version'] == "0.3":
            raise NotImplementedError

        gp_params['gp_orbit_prop_params'] = gp_orbit_prop_params
        gp_params['gp_orbit_link_params'] = gp_orbit_link_params
        gp_params['gp_general_params'] = gp_general_params
        gp_params['gp_instance_params'] = gp_instance_params
        gp_params['gp_data_rates_params'] = gp_data_rates_params
        gp_params['gp_other_params'] = gp_other_params
        gp_runner = GlobalPlannerRunner (gp_params)
        viz_outputs= gp_runner.run ()

        # define orbit prop outputs json
        output_json = {}
        output_json['version'] = OUTPUT_JSON_VER
        output_json['scenario_params'] = data['orbit_prop_inputs']['scenario_params']
        output_json['viz_data'] = viz_outputs
        output_json['update_time'] = datetime.utcnow().isoformat()

        return output_json


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description='orbit propagation')
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

    args = ap.parse_args()

    pr = PipelineRunner()

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
        "rs_s1_pickle": args.rs_s1_pickle,
        "rs_s2_pickle": args.rs_s2_pickle,
        "as_pickle": args.as_pickle,
        "file_params":  {'orbit_prop_inputs_file': args.prop_inputs_file.split('/')[-1]}
    }

    a = time.time()
    output = pr.run(data)
    b = time.time()
    with open('gp_outputs.json','w') as f:
        json.dump(output ,f)

    print('run time: %f'%(b-a))