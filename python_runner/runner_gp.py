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
        self.as_params = self.params['gp_general_params']['activity_scheduling_params']
        self.as_inst_params = self.params['gp_instance_params']['activity_scheduling_params']
        self.io_proc =GPProcessorIO(self.params)
        self.gp_plot =GPPlotting( self.params)

    def pickle_rtsel_stuff(self,routes_by_obs,all_stats,route_times_s,obs_indx,obs_winds,dlnk_winds_flat,ecl_winds,window_uid):

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
        pickle_name ='pickles/rs_%s_oi%d_%s' %( self.other_params['new_pickle_file_name_pre'],obs_indx,datetime.utcnow().isoformat().replace (':','_'))
        with open('%s.pkl' % ( pickle_name),'wb') as f:
            pickle.dump(  pickle_stuff,f)

    def pickle_actsc_stuff(self,routes_by_obs,ecl_winds,scheduled_routes,energy_usage,window_uid):

        pickle_stuff =  {}
        pickle_stuff['routes_by_obs'] = routes_by_obs
        pickle_stuff['ecl_winds'] = ecl_winds
        pickle_stuff['scheduled_routes'] = scheduled_routes
        pickle_stuff['energy_usage'] = energy_usage
        pickle_stuff['window_uid'] = window_uid
        pickle_name ='pickles/as_%s' %(datetime.utcnow().isoformat().replace (':','_'))
        with open('%s.pkl' % ( pickle_name),'wb') as f:
            pickle.dump(  pickle_stuff,f)

    def unpickle_rtsel_stuff( self):

        p = pickle.load (open ( self.pickle_params['route_selection_pickle'],'rb'))
        #  TODO:  uncommon this  if want to reload parameters from file
        # self.params = p['params']

        return p['routes_by_obs'],p['all_stats'],p['route_times_s'],p['obs_indx'],p['obs_winds'],p['dlnk_winds_flat'],p['ecl_winds'],p['window_uid']

    def unpickle_actsc_stuff( self):

        p = pickle.load (open ( self.pickle_params['act_scheduling_pickle'],'rb'))

        return p['routes_by_obs'],p['ecl_winds'],p['scheduled_routes'],p['energy_usage'],p['window_uid']

    

    def run_nominal_activity_scheduling( self, routes_by_obs,ecl_winds):
        gp_as = GPActivityScheduling ( self.params)

        # flatten the list of all routes, which currently has nested lists for each observation
        routes_flat = [rt for rts in routes_by_obs.values() for rt in rts]

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
        routes = gp_as.extract_utilized_routes ( copy_windows = True, verbose  = False)
        energy_usage = gp_as.extract_resource_usage(  decimation_factor =1)

        time_elapsed = t_b-t_a

        return  routes, energy_usage

    def calc_activity_scheduling_results ( self,obs_winds,dlnk_winds_flat,rs_routes_by_obs,sched_routes, energy_usage):
        gp_met = GPMetrics(self.params)

        total_collectible_DV_all_obs_winds = sum(obs.data_vol for winds in obs_winds for obs in winds)
        total_dlnkable_DV_all_dlnk_winds = sum(dlnk.data_vol for winds in dlnk_winds_flat for dlnk in winds)
        rs_output_routes = [rt for rts in rs_routes_by_obs.values() for rt in rts]
        total_throughput_DV_rs_routes = sum(sum(rt.data_vol for rt in rts) for obs, rts in rs_routes_by_obs.items())
        total_collectible_DV_rs_routes = sum(min(obs.data_vol,sum(rt.data_vol for rt in rts)) for obs, rts in rs_routes_by_obs.items())

        print('------------------------------')
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
        plot_outputs['obs_aoi_curves_by_targID'] = aoi_targ_stats['aoi_curves_by_targID_rs']
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

    def run_nominal_route_selection_v2( self,obs_winds,dlnk_winds_flat,xlnk_winds):
        gp_rs = gprsv2.GPDataRouteSelection ( self.params)

        print ('nominal route selection v2')

        obs_indx =0
        # dict of all routes, with obs as key
        routes_by_obs ={}
        all_stats =[]
        route_times_s =[]

        # this should be unique across all data routes
        dr_uid = 0

        for sat_indx in range( self.sat_params['num_sats']):
            for  index, obs in  enumerate ( obs_winds[sat_indx]):

                # obs.data_vol = 40000

                # if not (index == 0 and sat_indx == 2):
                #     continue

                print ("sat_indx")
                print (sat_indx)
                print ("obs")
                print ( index)

                # run the route selection algorithm
                t_a = time.time()
                routes,dr_uid = gp_rs.run_stage1(obs,dlnk_winds_flat,xlnk_winds, dr_uid, verbose = False)
                t_b = time.time()
                stats = gp_rs.get_stats(verbose=False )

                time_elapsed = t_b-t_a
                # end algorithm

                print ('len(routes)')
                print (len(routes))
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

                routes_by_obs[obs] = routes
                all_stats.append ( stats)
                route_times_s.append ( time_elapsed)

                obs_indx +=1

            #     if obs_indx >= 1:
            #         break

            # if obs_indx >= 1:
            #     break


        # gp_met = GPMetrics(self.params)
        # gp_met.assess_route_overlap( routes,verbose=False)

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
        for rts_indx, (obs,routes) in enumerate (routes_by_obs.items()):
            if rts_indx >= num_obs_to_plot:
                break

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
                    fig_name='plots/test_rs_o{0}_d{1}_6sat.pdf'.format (rts_indx,dlnk_indx)
                )


    def  plot_activity_scheduling_results ( self,routes_by_obs,routes,energy_usage,ecl_winds,metrics_plot_inputs):

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
        sats_to_include =  [sat_p['sat_id'] for sat_p in self.sat_orbit_params]
        # sats_to_include = [12,13,14,15,16]

        # sats_to_include =  range(20,30)
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
            plot_title = 'Scheduled Activities',
            plot_size_inches = (18,12),
            plot_include_dlnk_labels = self.as_params['plot_include_dlnk_labels'],
            plot_include_xlnk_labels = self.as_params['plot_include_xlnk_labels'],
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
            plot_title = 'Energy Utilization',
            plot_size_inches = (18,12),
            show=  False,
            fig_name='plots/test_energy.pdf'
        )

        targs_to_include = [targ['id'] for targ in self.obs_params['targets']]
        # targs_to_include = [0,3,4,7,8]
        targs_to_include = range(15)
        self.gp_plot.plot_obs_aoi(
            targs_to_include,
            metrics_plot_inputs['obs_aoi_curves_by_targID'],
            self.as_inst_params['start_utc_dt'],
            self.as_inst_params['end_utc_dt'],
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

            print('obs_winds')
            print(sum([len(p) for p in obs_winds]))
            print('dlnk_win')
            print(sum([len(p) for p in dlnk_winds]))
            print('xlnk_win')
            print(sum([len(xlnk_winds[i][j]) for i in  range( self.sat_params['num_sats']) for j in  range( self.sat_params['num_sats']) ]))

        #################################
        #  route selection stage
        #################################

        #  if  we are loading from file, do that
        if self.pickle_params['unpickle_pre_route_selection']:
            routes_by_obs,all_stats,route_times_s,obs_indx,obs_winds,dlnk_winds_flat,ecl_winds,window_uid = self.unpickle_rtsel_stuff()

        #  otherwise run route selection
        else:
            routes_by_obs,all_stats,route_times_s, obs_indx, dr_uid  =  self.run_nominal_route_selection_v2(obs_winds,dlnk_winds_flat,xlnk_winds)
            # routes_by_obs,all_stats,route_times_s, obs_indx, weights_tups  =  self.run_test_route_selection(obs_winds,dlnk_winds_flat,xlnk_winds)

        if self.pickle_params['pickle_route_selection_results']:
            self.pickle_rtsel_stuff(routes_by_obs,all_stats,route_times_s,obs_indx,obs_winds,dlnk_winds_flat,ecl_winds,window_uid)

        # run stage 2. todo:  move this elsewhere
        all_routes = []

        # import ipdb
        # ipdb.set_trace()

        for obs,routes in routes_by_obs.items():
            all_routes += routes

        print('len(all_routes)')
        print(len(all_routes))

        gp_met = GPMetrics(self.params)
        # t_a = time.time()
         # need to figure out how to window this or something so that we don't have to compare to every other data route - that's horribly expensive
        gp_met.assess_route_overlap( all_routes,verbose=True)
        # t_b = time.time()
        # time_elapsed = t_b-t_a
        # print('time_elapsed')
        # print(time_elapsed)
        
        #################################
        # route selection output stage
        #################################
        
        print('route selection output stage')

        if self.rs_general_params['plot_route_selection_results']:
            self.plot_route_selection_results (routes_by_obs,dlnk_winds_flat,xlnk_winds_flat,num_obs_to_plot = 5)

        #################################
        #  activity scheduling stage
        #################################

        if not self.as_params['run_activity_scheduling']:
            return None

        #  explicitly validate routes 
        # for routes in routes_by_obs:
        #     for dr in routes:
        #         dr.validate_route()

        #  if  we are loading from file, do that
        if self.pickle_params['unpickle_pre_act_scheduling']:
            routes_by_obs,ecl_winds,scheduled_routes,energy_usage,window_uid = self.unpickle_actsc_stuff()
        else:
            scheduled_routes,energy_usage = self.run_nominal_activity_scheduling(routes_by_obs,ecl_winds)

        # if we are saving to file, do that
        if self.pickle_params['pickle_act_scheduling_results']:
            self.pickle_actsc_stuff(routes_by_obs,ecl_winds,scheduled_routes,energy_usage,window_uid)

        metrics_plot_inputs = self.calc_activity_scheduling_results (obs_winds,dlnk_winds_flat,routes_by_obs,scheduled_routes, energy_usage)

        #################################
        #  Activity scheduling output stage
        #################################
        
        print('activity scheduling output stage')

        if self.as_params['plot_activity_scheduling_results']:
            self.plot_activity_scheduling_results(routes_by_obs,scheduled_routes,energy_usage,ecl_winds,metrics_plot_inputs)
          

        outputs = None
        # outputs= self.io_proc.make_sat_history_outputs (sel_obs_winds_flat, sel_xlnk_winds_flat, sel_dlnk_winds_flat, link_info_by_wind)

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

        if orbit_prop_inputs['version'] == "0.3": 
            # do some useful transformations while preserving the structure of the inputs ( important for avoiding namespace clashes)
            orbit_prop_inputs['scenario_params']['start_utc_dt'] = tt.iso_string_to_dt ( orbit_prop_inputs['scenario_params']['start_utc'])
            orbit_prop_inputs['scenario_params']['end_utc_dt'] = tt.iso_string_to_dt ( orbit_prop_inputs['scenario_params']['end_utc'])
            orbit_prop_inputs['sat_params']['num_sats'] = orbit_prop_inputs['sat_params']['num_satellites']
            orbit_prop_inputs['gs_params']['num_gs'] = orbit_prop_inputs['gs_params']['num_stations']
            orbit_prop_inputs['sat_params']['pl_data_rate'] = orbit_prop_inputs['sat_params']['payload_data_rate_Mbps']
            # orbit_prop_inputs['sat_orbit_params'], dummy = io_tools.unpack_sat_entry_list( orbit_prop_inputs['sat_orbit_params'],force_duplicate =  True)
            orbit_prop_inputs['sat_params']['power_params'], all_sat_ids1 = io_tools.unpack_sat_entry_list( orbit_prop_inputs['sat_params']['power_params'])
            orbit_prop_inputs['sat_params']['initial_state'], all_sat_ids2 = io_tools.unpack_sat_entry_list( orbit_prop_inputs['sat_params']['initial_state'])

            #  check if  we saw the same list of satellite IDs from each unpacking. if not that's a red flag that the inputs could be wrongly specified
            if all_sat_ids1 != all_sat_ids2:
                raise Exception('Saw differing sat ID orders')

            #  grab the list for satellite ID order.  if it's "default", we will create it and save it for future use here
            sat_id_order=orbit_prop_inputs['sat_params']['sat_id_order']
            #  make the satellite ID order. if the input ID order is default, then will assume that the order is the same as all of the IDs found in the power parameters
            sat_id_order = io_tools.make_and_validate_sat_id_order(sat_id_order,orbit_prop_inputs['sat_params']['num_sats'],all_sat_ids1)
            orbit_prop_inputs['sat_params']['sat_id_order'] = sat_id_order

            orbit_prop_inputs['sat_params']['power_params_sorted'] = io_tools.sort_input_params_by_sat_indcs(orbit_prop_inputs['sat_params']['power_params'],sat_id_order)
            orbit_prop_inputs['sat_params']['initial_state_sorted'] = io_tools.sort_input_params_by_sat_indcs(orbit_prop_inputs['sat_params']['initial_state'],sat_id_order)
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
        if not data_rates_inputs['version'] == "0.2": 
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

    args = ap.parse_args()

    pr = PipelineRunner()

    # with open(os.path.join(REPO_BASE,'crux/config/examples/orbit_prop_inputs_ex.json'),'r') as f:
    with open(os.path.join(REPO_BASE, args.prop_inputs_file),'r') as f:
        orbit_prop_inputs = json.load(f)

    with open(os.path.join(REPO_BASE,args.data_rates_file),'r') as f:
        data_rates_inputs = json.load(f)
        
    with open(os.path.join(REPO_BASE, args.link_inputs_file),'r') as f:
        orbit_link_inputs = json.load(f)

    with open(os.path.join(REPO_BASE,'crux/config/examples/gp_general_params_inputs_ex.json'),'r') as f:
        gp_general_params_inputs = json.load(f)

    with open(os.path.join(REPO_BASE,'crux/config/examples/gp_instance_params_inputs_ex.json'),'r') as f:
        gp_instance_params_inputs = json.load(f)

    data = {
        # "orbit_prop_data": orbit_prop_data,
        "orbit_prop_inputs": orbit_prop_inputs,
        "orbit_link_inputs": orbit_link_inputs,
        "gp_general_params_inputs": gp_general_params_inputs,
        "gp_instance_params_inputs": gp_instance_params_inputs,
        # "viz_params": viz_params,
        "data_rates_inputs": data_rates_inputs,
        "file_params":  {'orbit_prop_inputs_file': args.prop_inputs_file.split('/')[-1]}
    }

    a = time.time()
    output = pr.run(data)
    b = time.time()
    with open('gp_outputs.json','w') as f:
        json.dump(output ,f)

    print('run time: %f'%(b-a))