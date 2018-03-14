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


#  local repo includes. todo:  make this less hackey
sys.path.append ('..')
from circinus_tools  import time_tools as tt
from circinus_tools  import io_tools
from gp_tools.input_processing import GPProcessorIO
from gp_tools.gp_plotting import GPPlotting
from gp_tools.gp_route_selection import GPDataRouteSelection
from gp_tools.gp_activity_scheduling import GPActivityScheduling

# TODO: remove this line if not needed
from gp_tools.custom_activity_window import ObsWindow

REPO_BASE = os.path.abspath(os.pardir)  # os.pardir aka '..'

OUTPUT_JSON_VER = '0.1'

class GlobalPlannerRunner:
    """easy interface for running the global planner scheduling algorithm"""

    def __init__(self, params):
        self.params = params
        self.general_params = params['general_params']
        self.pickle_params = params['pickle_params']
        self.other_params = params['other_params']
        self.route_selection_params = params['route_selection_params']
        self.activity_scheduling_params = params['activity_scheduling_params']
        self.plot_params = params['plot_params']
        self.io_proc =GPProcessorIO( dict(list (self.general_params.items()) +  list (self.other_params.items()) ))
        self.gp_plot =GPPlotting( self.plot_params)

    def pickle_stuff(self,all_routes,all_routes_obs,all_stats,route_times_s,obs_indx):

        pickle_stuff =  {}
        pickle_stuff['all_routes'] = all_routes
        pickle_stuff['all_routes_obs'] = all_routes_obs
        pickle_stuff['all_stats'] = all_stats
        pickle_stuff['route_times_s'] = route_times_s
        pickle_stuff['params'] =  self.params
        pickle_stuff['obs_indx'] = obs_indx
        pickle_name ='%s_obsindx%d_%s' %( self.general_params['new_pickle_file_name_pre'],obs_indx,datetime.utcnow().isoformat().replace (':','_'))
        with open('%s.pkl' % ( pickle_name),'wb') as f:
            pickle.dump(  pickle_stuff,f)

    def unpickle_stuff( self):

        p = pickle.load (open ( self.pickle_params['route_selection_pickle'],'rb'))
        self.params = p['params']

        return p['all_routes'],p['all_routes_obs'],p['all_stats'],p['route_times_s'],p['obs_indx']

    

    def run_nominal_activity_scheduling( self, all_routes):
        
        gp_as = GPActivityScheduling ( dict(list (self.general_params.items()) + list (self.activity_scheduling_params. items())))

        # flatten the list of all routes, which currently has nested lists for each observation
        routes_flat = [item for sublist in all_routes for item in sublist]

        print('make activity scheduling model')
        gp_as.make_model (routes_flat, verbose = True)
        stats =gp_as.get_stats (verbose = True)
        print('solve activity scheduling')
        t_a = time.time()
        gp_as.solve ()
        t_b = time.time()
        # gp_as.print_sol ()
        print('extract_routes')
        routes = gp_as. extract_utilized_routes ( verbose  = False)

        print('len(routes)')
        print(len(routes))

        time_elapsed = t_b-t_a

        return  routes

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

        gp_ps.make_model (obs,dlnk_winds_flat,xlnk_winds, obj_weights, verbose = True)
        stats =gp_ps.get_stats (verbose = True)
        t_a = time.time()
        gp_ps.solve ()
        t_b = time.time()
        gp_ps.print_sol ()
        routes,dr_uid = gp_ps. extract_routes ( dr_uid,verbose  = True)

        time_elapsed = t_b-t_a

        return routes,obs,stats,time_elapsed,dr_uid

    def run_nominal_route_selection( self,obs_winds,dlnk_winds_flat,xlnk_winds):
        gp_ps = GPDataRouteSelection ( dict(list (self.general_params.items()) + list (self.route_selection_params. items())))

        print ('nominal route selection')

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

        for sat_indx in range( self.general_params['num_sats']):
            for  index, obs in  enumerate ( obs_winds[sat_indx]):

                print ("sat_indx")
                print (sat_indx)
                print ("obs")
                print ( index)

                routes,obs,stats,time_elapsed,dr_uid = self.run_route_selection(gp_ps,obs,dlnk_winds_flat,xlnk_winds,obj_weights,dr_uid)

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

            #     if obs_indx >= 2:
            #         break

            # if obs_indx >= 2:
            #     break

        return all_routes,all_routes_obs,all_stats,route_times_s,obs_indx

    def  setup_test( self,obs_winds,dlnk_winds_flat,xlnk_winds):
        obs_winds_sel = [[] for sat_indx in range (self.general_params['num_sats'])]
        # obs_winds_sel[13].append ( obs_winds[13][0])
        # obs_winds_sel[19].append ( obs_winds[19][0])
        # obs_winds_sel[19].append ( obs_winds[19][1])
        # obs_winds_sel[20].append ( obs_winds[20][6])
        # obs_winds_sel[26].append ( obs_winds[26][0])

        if False:
            self.gp_plot.plot_winds(
                range (self.general_params['num_sats']),
                # [13,19,20,26],
                obs_winds_sel,
                [],
                dlnk_winds_flat, 
                [],
                [],
                {},
                self.general_params['start_utc_dt'],
                self.general_params['end_utc_dt'],
                plot_title = 'all obs, dlnks', 
                plot_size_inches = (18,12),
                plot_include_labels = True,
                show=  False,
                fig_name='plots/temp1.pdf'
            )

        return obs_winds_sel

        
    def run_test_route_selection( self,obs_winds,dlnk_winds_flat,xlnk_winds):
        gp_ps = GPDataRouteSelection ( dict(list (self.general_params.items()) + list (self.route_selection_params. items())))

        total_dv_weights = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] 
        num_paths_sel_weights = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1] 
        latency_sf_weights = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0] 
        weights_tups = zip(total_dv_weights,num_paths_sel_weights,latency_sf_weights)

        print ('test route selection')

        # obs_winds_sel =  self.setup_test(obs_winds,dlnk_winds_flat,xlnk_winds)
        obs_winds_sel = [[] for sat_indx in range (self.general_params['num_sats'])]
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

                routes,obs,stats,time_elapsed = self.run_route_selection(gp_ps,obs,dlnk_winds_flat,xlnk_winds,obj_weights)

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

    def plot_route_selection_results ( self,obs_routes,dlnk_winds_flat,xlnk_winds_flat):

        for rts_indx, routes in enumerate (obs_routes):

            # TODO:  this stuff needs to be  changed
            sel_obs_winds_flat, sel_dlnk_winds_flat, \
            sel_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind = self.io_proc.extract_flat_windows (routes)

            obs = None
            for sat_indx in  range (self.general_params['num_sats']):
                for obs in sel_obs_winds_flat[sat_indx]:
                    if obs:
                        break

            sats_to_include =  range (self.general_params['num_sats'])
            # sats_to_include = [0,9,21,22,23]

            #  plot the selected down links and cross-links
            self.gp_plot.plot_winds(
                sats_to_include,
                sel_obs_winds_flat,
                dlnk_winds_flat,
                sel_dlnk_winds_flat, 
                xlnk_winds_flat,
                sel_xlnk_winds_flat,
                route_indcs_by_wind,
                obs.start,
                obs.start + timedelta( seconds= self.route_selection_params['wind_filter_duration_s']),
                # self.general_params['start_utc_dt'],
                # self.general_params['start_utc_dt'] + timedelta( seconds= self.route_selection_params['wind_filter_duration_s']),
                # self.general_params['end_utc_dt']-timedelta(minutes=200),
                plot_title = 'Route Plot', 
                plot_size_inches = (18,12),
                plot_include_labels = self.route_selection_params['plot_include_labels'],
                show= False,
                fig_name='plots/obs_winds_20_6_rts{0}.pdf'.format (rts_indx)
            )


    def  plot_activity_scheduling_results ( self,obs_routes,routes):

        # do a bunch of stuff to extract the windows from all of the routes as indexed by observation
        # note that this stuff is not thewindows from the scheduled routes, but rather the windows from all the route selected in route selection
        #  start
        sel_obs_winds_flat = [set() for  sat_indx  in range  (self.general_params['num_sats'])]
        sel_dlnk_winds_flat = [set() for sat_indx  in range (self.general_params['num_sats'])]
        sel_xlnk_winds_flat = [set() for sat_indx  in range (self.general_params['num_sats'])]        

        for rts_indx, rts in enumerate (obs_routes):
            obs_winds_rt, dlnk_winds_rt, \
            xlnk_winds_rt, _, _ = self.io_proc.extract_flat_windows (rts)

            for sat_indx in range  (self.general_params['num_sats']):
                [sel_obs_winds_flat[sat_indx] .add( wind)  for wind in obs_winds_rt[sat_indx]]
                [sel_dlnk_winds_flat[sat_indx] .add(wind ) for wind in dlnk_winds_rt[sat_indx]]
                [sel_xlnk_winds_flat[sat_indx] .add(wind ) for wind in xlnk_winds_rt[sat_indx]]

        for sat_indx in range  (self.general_params['num_sats']):
            sel_obs_winds_flat[sat_indx] = list(sel_obs_winds_flat[sat_indx])
            sel_dlnk_winds_flat[sat_indx] = list(sel_dlnk_winds_flat[sat_indx])
            sel_xlnk_winds_flat[sat_indx] = list(sel_xlnk_winds_flat[sat_indx])
            sel_obs_winds_flat[sat_indx].sort(key=lambda x: x.start)
            sel_dlnk_winds_flat[sat_indx].sort(key=lambda x: x.start)
            sel_xlnk_winds_flat[sat_indx].sort(key=lambda x: x.start)
        # end

        sched_obs_winds_flat, sched_dlnk_winds_flat, \
        sched_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind = self.io_proc.extract_flat_windows (routes,copy_windows=True)

        #  update the window beginning and end times based upon their amount of scheduled data volume
        for sat_indx in range(self.general_params['num_sats']):
            for wind in sched_obs_winds_flat[sat_indx] + sched_dlnk_winds_flat[sat_indx] + sched_xlnk_winds_flat[sat_indx]:
                print(wind)
                print(wind.data_vol)
                print(wind.scheduled_data_vol)
                if not  type (wind) == ObsWindow:
                    print(link_info_by_wind[wind])
                    print(route_indcs_by_wind[wind])
                wind.update_duration_from_scheduled_dv ()

        sats_to_include =  range (self.general_params['num_sats'])
        # sats_to_include = [0,9,21,22,23]


        #  plot the selected down links and cross-links
        self.gp_plot.plot_winds(
            sats_to_include,
            sel_obs_winds_flat,
            sched_obs_winds_flat,
            sel_dlnk_winds_flat,
            sched_dlnk_winds_flat, 
            sel_xlnk_winds_flat,
            sched_xlnk_winds_flat,
            route_indcs_by_wind,
            self.general_params['start_utc_dt'],
            # self.general_params['start_utc_dt'] + timedelta( seconds= self.route_selection_params['wind_filter_duration_s']),
            self.general_params['end_utc_dt'],
            plot_title = 'Scheduled Activities',
            plot_size_inches = (18,12),
            plot_include_labels = self.activity_scheduling_params['plot_include_labels'],
            show=  False,
            fig_name='plots/test.pdf'
        )

    def validate_unique_windows( self,obs_winds,dlnk_winds_flat,xlnk_winds):
        all_wind_ids = []

        for indx in range(len(obs_winds)):
            for wind in obs_winds[indx]:
                if wind.window_ID in all_wind_ids:
                    raise Exception('Found a duplicate unique window ID where it should not have been possible')
                all_wind_ids.append(wind.window_ID)

        for indx in range(len(dlnk_winds_flat)):
            for wind in dlnk_winds_flat[indx]:
                if wind.window_ID in all_wind_ids:
                    raise Exception('Found a duplicate unique window ID where it should not have been possible')
                all_wind_ids.append(wind.window_ID)

        for indx in range(len(xlnk_winds)):
            for indx_2 in range(len(xlnk_winds[indx])):
                for wind in xlnk_winds[indx][indx_2]:
                    if wind.window_ID in all_wind_ids:
                        raise Exception('Found a duplicate unique window ID where it should not have been possible')
                    all_wind_ids.append(wind.window_ID)

    def run( self):

        #################################
        #  parse inputs, if desired
        #################################

        if not self.pickle_params['unpickle_pre_route_selection']:
            # parse the inputs into activity windows
            window_uid = 0
            obs_winds, window_uid =self.io_proc.import_obs_winds(window_uid)
            dlnk_winds, dlnk_winds_flat, window_uid =self.io_proc.import_dlnk_winds(window_uid)
            xlnk_winds, xlnk_winds_flat, window_uid =self.io_proc.import_xlnk_winds(window_uid)

            # important to check this because window unique IDs are used as hashes in dictionaries in the scheduling code
            self.validate_unique_windows(obs_winds,dlnk_winds_flat,xlnk_winds)

            # todo:  probably ought to delete the input times and rates matrices to free up space

            print('obs_winds')
            print(sum([len(p) for p in obs_winds]))
            print('dlnk_win')
            print(sum([len(p) for p in dlnk_winds]))
            print('xlnk_win')
            print(sum([len(xlnk_winds[i][j]) for i in  range( self.general_params['num_sats']) for j in  range( self.general_params['num_sats']) ]))

        #################################
        #  route selection stage
        #################################

        #  if  we are loading from file, do that
        if self.pickle_params['unpickle_pre_route_selection']:
            obs_routes,all_routes_obs,all_stats,route_times_s,obs_indx = self.unpickle_stuff()

        #  otherwise run route selection
        else:
            obs_routes,all_routes_obs,all_stats,route_times_s, obs_indx  =  self.run_nominal_route_selection(obs_winds,dlnk_winds_flat,xlnk_winds)
            # obs_routes,all_routes_obs,all_stats,route_times_s, obs_indx, weights_tups  =  self.run_test_route_selection(obs_winds,dlnk_winds_flat,xlnk_winds)

        if self.pickle_params['pickle_route_selection_results']:
            self.pickle_stuff(obs_routes,all_routes_obs,all_stats,route_times_s,obs_indx)
        
        #################################
        #  activity scheduling stage
        #################################

        #  explicitly validate routes 
        # for routes in obs_routes:
        #     for dr in routes:
        #         dr.validate_route()


        scheduled_routes = self.run_nominal_activity_scheduling(obs_routes)

        #################################
        #   output stage
        #################################
        
        print('output stage')

        if self.route_selection_params['plot_route_selection_results']:
            self.plot_route_selection_results (obs_routes,dlnk_winds_flat,xlnk_winds_flat)

        if self.activity_scheduling_params['plot_activity_scheduling_results']:
            self.plot_activity_scheduling_results(obs_routes,scheduled_routes)
          

        outputs = None
        # outputs= self.io_proc.make_sat_history_outputs (sel_obs_winds_flat, sel_xlnk_winds_flat, sel_dlnk_winds_flat, link_info_by_wind)

        return outputs


class PipelineRunner:

    def run(self, data):
        """

        """

        orbit_prop_inputs = data['orbit_prop_inputs']
        gp_params_inputs = data['gp_params_inputs']
        data_rates_output = data['data_rates_output']
        file_params = data.get ('file_params',{})

        gp_params = {}
        gp_general_params = {}
        gp_general_params['new_pickle_file_name_pre']  = 'pickles/'  + file_params.get ('orbit_prop_inputs_file' ,'default.json').split ('.')[0]

        if orbit_prop_inputs['version'] == "0.3": 
            
            gp_general_params['start_utc_dt'] = tt.iso_string_to_dt ( orbit_prop_inputs['scenario_params']['start_utc'])
            gp_general_params['end_utc_dt'] = tt.iso_string_to_dt ( orbit_prop_inputs['scenario_params']['end_utc'])
            gp_general_params['timestep_s'] = orbit_prop_inputs['scenario_params']['timestep_s']
            gp_general_params['num_sats'] = orbit_prop_inputs['sat_params']['num_satellites']
            gp_general_params['num_gs'] = orbit_prop_inputs['gs_params']['num_stations']
            gp_general_params['num_targets'] = orbit_prop_inputs['obs_params']['num_targets']
            gp_general_params['pl_data_rate'] = orbit_prop_inputs['sat_params']['payload_data_rate_Mbps']
            gp_general_params['power_params'], all_sat_ids1 = io_tools.unpack_sat_entry_list( orbit_prop_inputs['sat_params']['power_params'])
            gp_general_params['initial_state'], all_sat_ids2 = io_tools.unpack_sat_entry_list( orbit_prop_inputs['sat_params']['initial_state'])
            #  check if  we saw the same list of satellite IDs from each unpacking. if not that's a red flag that the inputs could be wrongly specified
            if all_sat_ids1 != all_sat_ids2:
                raise Exception('Saw differing sat ID orders')

            #  grab the list for satellite ID order.  if it's "default", we will create it and save it for future use here
            sat_id_order=orbit_prop_inputs['sat_params']['sat_id_order']
            sat_id_order = io_tools.make_and_validate_sat_id_order(sat_id_order,gp_general_params['num_sats'],all_sat_ids1)
            gp_general_params['sat_id_order'] = sat_id_order


        if gp_params_inputs['version'] == "0.1": 
            gp_params['other_params'] = gp_params_inputs['other_params']
            gp_params['route_selection_params'] = gp_params_inputs['route_selection_params']
            gp_params['activity_scheduling_params'] = gp_params_inputs['activity_scheduling_params']
            gp_params['pickle_params'] = gp_params_inputs['pickle_params']
            gp_params['plot_params'] = gp_params_inputs['plot_params']

        if data_rates_output['version'] == "0.1": 
            gp_general_params['obs_times'] = data_rates_output['accesses_data_rates']['obs_times']
            gp_general_params['dlnk_times'] = data_rates_output['accesses_data_rates']['dlnk_times']
            gp_general_params['dlnk_rates'] = data_rates_output['accesses_data_rates']['dlnk_rates']
            gp_general_params['xlnk_times'] = data_rates_output['accesses_data_rates']['xlnk_times']
            gp_general_params['xlnk_rates'] = data_rates_output['accesses_data_rates']['xlnk_rates']

        gp_params['general_params'] = gp_general_params
        gp_runner = GlobalPlannerRunner (gp_params)
        viz_outputs= gp_runner.run ()

        # define orbit prop outputs json
        output_json = {}
        output_json['version'] = OUTPUT_JSON_VER
        output_json['scenario_params'] = orbit_prop_inputs['scenario_params']
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

    args = ap.parse_args()

    pr = PipelineRunner()

    # with open(os.path.join(REPO_BASE,'crux/config/examples/orbit_prop_inputs_ex.json'),'r') as f:
    with open(os.path.join(REPO_BASE, args.prop_inputs_file),'r') as f:
        orbit_prop_inputs = json.load(f)

    with open(os.path.join(REPO_BASE,args.data_rates_file),'r') as f:
        data_rates_output = json.load(f)

    with open(os.path.join(REPO_BASE,'crux/config/examples/gp_params_inputs_ex.json'),'r') as f:
        gp_params_inputs = json.load(f)

    data = {
        # "orbit_prop_data": orbit_prop_data,
        "orbit_prop_inputs": orbit_prop_inputs,
        "gp_params_inputs": gp_params_inputs,
        # "viz_params": viz_params,
        "data_rates_output": data_rates_output,
        "file_params":  {'orbit_prop_inputs_file': args.prop_inputs_file.split('/')[-1]}
    }

    a = time.time()
    output = pr.run(data)
    b = time.time()
    with open('gp_outputs.json','w') as f:
        json.dump(output ,f)

    print('run time: %f'%(b-a))