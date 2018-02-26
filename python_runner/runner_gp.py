#! /usr/bin/env python

##
# Python runner for orbit visualization pipeline
# @author Kit Kennedy
# 

import time
import os.path
import json
from datetime import datetime
import sys
import argparse

#  local repo includes. todo:  make this less hackey
sys.path.append ('..')
from circinus_tools  import time_tools as tt
from gp_tools.input_processing import GPProcessorIO
from gp_tools.gp_route_selection import GPDataRouteSelection


REPO_BASE = os.path.abspath(os.pardir)  # os.pardir aka '..'

OUTPUT_JSON_VER = '0.1'

class GlobalPlannerRunner:
    """easy interface for running the global planner scheduling algorithm"""

    def __init__(self, params):
        self.params = params
        self.io_proc =GPProcessorIO( params)

    def run( self):

        # parse the inputs into activity windows
        obs_winds, obs_window_id =self.io_proc.import_obs_winds()
        dlink_winds, dlink_winds_flat, dlnk_window_id =self.io_proc.import_dlnk_winds()
        xlink_winds, xlink_winds_flat, xlnk_window_id =self.io_proc.import_xlnk_winds()

        # for j in dlink_winds:
        #     for i in j:
        #         i.print_self ()

        # todo:  probably ought to delete the input times and rates matrices to free up space

        print('obs_winds')
        print(sum([len(p) for p in obs_winds]))
        print('dlink_win')
        print(sum([len(p) for p in dlink_winds]))
        print('xlink_win')
        print(sum([len(xlink_winds[i][j]) for i in  range( self.params['num_sats']) for j in  range( self.params['num_sats']) ]))

        gp_ps = GPDataRouteSelection ( self.params)
        gp_ps.make_model (obs_winds[5][0],dlink_winds_flat,xlink_winds)
        gp_ps.solve ()
        print ('obs_winds[5][0].sat_indx')
        print (obs_winds[5][0].sat_indx)
        print ('obs_winds[5][0].data_vol')
        print (obs_winds[5][0].data_vol)
        print ('obs_winds[5][0].duration')
        print (obs_winds[5][0].duration)
        print ("(obs_winds[5][0].end-self.params['start_utc_dt']).total_seconds ()")  
        print ( (obs_winds[5][0].end-self.params['start_utc_dt']).total_seconds ())  
        gp_ps.print_sol ()
        routes = gp_ps. extract_routes ( verbose  = True)


        sel_obs_winds_flat, sel_xlnk_winds_flat, \
        sel_dlnk_winds_flat, link_info_by_wind = self.io_proc.extract_flat_windows ( routes)

        outputs= self.io_proc.make_sat_history_outputs (sel_obs_winds_flat, sel_xlnk_winds_flat, sel_dlnk_winds_flat, link_info_by_wind)

        return outputs


class PipelineRunner:

    def run(self, data):
        """

        """

        orbit_prop_inputs = data['orbit_prop_inputs']
        gp_params_inputs = data['gp_params_inputs']
        data_rates_output = data['data_rates_output']

        gp_params = {}

        if orbit_prop_inputs['version'] == "0.3": 
            
            gp_params['start_utc_dt'] = tt.iso_string_to_dt ( orbit_prop_inputs['scenario_params']['start_utc'])
            gp_params['end_utc_dt'] = tt.iso_string_to_dt ( orbit_prop_inputs['scenario_params']['end_utc'])
            gp_params['timestep_s'] = orbit_prop_inputs['scenario_params']['timestep_s']
            gp_params['num_sats'] = orbit_prop_inputs['sat_params']['num_satellites']
            gp_params['num_gs'] = orbit_prop_inputs['gs_params']['num_stations']
            gp_params['pl_data_rate'] = orbit_prop_inputs['sat_params']['payload_data_rate_Mbps']

        if gp_params_inputs['version'] == "0.1": 
            gp_params['targ_ignore_list'] = gp_params_inputs['targ_ignore_list']
            gp_params['gs_ignore_list'] = gp_params_inputs['gs_ignore_list']
            gp_params['min_allowed_dv_dlnk'] = gp_params_inputs['min_allowed_dv_dlnk_Mb']
            gp_params['min_allowed_dv_xlnk'] = gp_params_inputs['min_allowed_dv_xlnk_Mb']
            gp_params['route_selection_num_paths'] = gp_params_inputs['route_selection_num_paths']
            gp_params['route_selection_min_path_dv'] = gp_params_inputs['route_selection_min_path_dv_Mb']
            gp_params['solver_max_runtime'] = gp_params_inputs['solver_max_runtime_s']

        if data_rates_output['version'] == "0.1": 
            gp_params['obs_times'] = data_rates_output['accesses_data_rates']['obs_times']
            gp_params['dlnk_times'] = data_rates_output['accesses_data_rates']['dlnk_times']
            gp_params['dlnk_rates'] = data_rates_output['accesses_data_rates']['dlnk_rates']
            gp_params['xlnk_times'] = data_rates_output['accesses_data_rates']['xlnk_times']
            gp_params['xlnk_rates'] = data_rates_output['accesses_data_rates']['xlnk_rates']

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
        "data_rates_output": data_rates_output
    }

    a = time.time()
    output = pr.run(data)
    b = time.time()
    with open('gp_outputs.json','w') as f:
        json.dump(output ,f)

    print('run time: %f'%(b-a))