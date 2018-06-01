# contains model and Solver for global planner activity scheduling capability
# 
# @author Kit Kennedy
#

from  datetime import timedelta
from copy import  deepcopy
from math import ceil

from circinus_tools  import time_tools as tt
from circinus_tools  import  constants as const
from circinus_tools  import io_tools
from circinus_tools.scheduling.custom_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow
from circinus_tools.scheduling.schedule_objects import Dancecard
from circinus_tools.scheduling.routing_objects import DataMultiRoute
from circinus_tools.sat_state_tools import propagate_sat_ES
from circinus_tools.scheduling.formulation.agent_scheduler import AgentScheduling

class GPActivityScheduling(AgentScheduling):
    """Superclass for GP activity scheduling"""
    
    def __init__(self,gp_params):
        """initializes based on parameters
        
        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """
        super().__init__()

        scenario_params = gp_params['orbit_prop_params']['scenario_params']
        sat_params = gp_params['orbit_prop_params']['sat_params']
        as_params = gp_params['gp_general_params']['activity_scheduling_params']
        gp_inst_planning_params = gp_params['gp_instance_params']['planning_params']

        self.gp_inst_planning_params = gp_inst_planning_params
        self.gp_agent_ID = gp_params['gp_instance_params']['gp_agent_ID']

        self.scenario_timestep_s = scenario_params['timestep_s']
        
        self.latency_calculation_params = gp_params['gp_general_params']['other_params']['latency_calculation']
        self.solver_name =as_params['solver_name']
        self.solver_params =as_params['solver_params']
        self.num_sats=sat_params['num_sats']
        self.resource_delta_t_s  =as_params['resource_delta_t_s']
        self.enforce_energy_storage_constr  =as_params['enforce_energy_storage_constr']
        self.enforce_data_storage_constr  =as_params['enforce_data_storage_constr']

        # this is the minimum obs dv that must be downlinked for an obs (data route/obs dlnk) in order for it to count it towards objective terms (other than total dv)
        self.min_obs_dv_dlnk_req =as_params['min_obs_dv_dlnk_req_Mb']
        #  the amount of extra utilization a fixed route is allowed in the model.
        #  the "effectively zero" number.
        self.dv_epsilon = as_params['dv_epsilon_Mb']
        self.fixed_utilization_epsilon = as_params['fixed_utilization_epsilon']

        # See docs for this in agent_scheduler superclass
        self.min_latency_for_sf_1_mins =as_params['min_latency_for_sf_1_mins']

        # notes on the planning window:
        # -  new data routes are filtered such that all of their activity windows must fall completely within planning_fixed_end_dt and planning_end_dt. 
        # - existing data routes are filtered such that any of their activity windows must fall partially within planning_start_dt and planning_end_dt. 
        # - the wider filter for existing data routes means that we capture all activity windows that are already scheduled and exert the constraints that are present for them. existing data routes may only be "scheduled downwards"; that is, their utilization decreases. this means that we only have to care about the activity windows that are within the planning window for these routes, because the the system is as or less constrained when utilization is lowered ( i.e., we don't pretend that data volume came out of thin air in the past in order to increase the utilization above what was already scheduled for an existing data route)
        self.planning_start_dt  = gp_inst_planning_params['planning_start_dt']
        # planning_fixed_end is the time before which any pre-existing data routes will be considered "already scheduled for good", and their utilization will be fixed, unchangeable. This reckoning is based on the beginning of the route (start of the obs window) - if the start of the obs window falls before planning_fixed_end, then the route will be considered fixed. Also, no activity windows from new routes may be scheduled before this time.
        self.planning_fixed_end_dt  = gp_inst_planning_params['planning_fixed_end_dt']
        self.planning_end_obs_dt = gp_inst_planning_params['planning_end_obs_dt']
        self.planning_end_xlnk_dt = gp_inst_planning_params['planning_end_xlnk_dt']
        self.planning_end_dlnk_dt  = gp_inst_planning_params['planning_end_dlnk_dt']
        self.planning_end_dt  = self.planning_end_dlnk_dt

        self.resource_margin_obj_num_timepoints = as_params['resource_margin_obj_num_timepoints']

        self.obj_weights =as_params['obj_weights']

        self.use_symmetric_xlnk_windows = gp_params['gp_general_params']['other_params']['use_symmetric_xlnk_windows']

        self.power_params = sat_params['power_params_sorted']
        self.data_storage_params = sat_params['data_storage_params_sorted']
        self.sats_state = gp_params['gp_instance_params']['sats_state_sorted']
        self.sat_activity_params = sat_params['activity_params']


        self.sats_dmin_Mb = [1000*ds_params['data_storage_Gbit']['d_min'][ds_params['storage_option']] for ds_params in self.data_storage_params]
        self.sats_dmax_Mb = [1000*ds_params['data_storage_Gbit']['d_max'][ds_params['storage_option']] for ds_params in self.data_storage_params]

        self.energy_unit = "Wh"  # watt hours

        # these lists are in order of satellite index because we've sorted 
        self.sats_init_estate_Wh = [sat_state['batt_e_Wh'] for sat_state in self.sats_state]
        self.sats_edot_by_mode_W = []
        self.sats_emin_Wh = []
        self.sats_emax_Wh = []
        for p_params in self.power_params:
            sat_edot_by_mode,sat_batt_storage,power_units = io_tools.parse_power_consumption_params(p_params)

            if not power_units['power_consumption'] == 'W':
                raise NotImplementedError
            if not power_units['battery_storage'] == 'Wh':
                raise NotImplementedError

            self.sats_edot_by_mode_W.append (sat_edot_by_mode)
            self.sats_emin_Wh.append (sat_batt_storage['e_min'])
            self.sats_emax_Wh.append (sat_batt_storage['e_max'])

        self.min_act_duration_s = {
            ObsWindow: self.sat_activity_params['min_duration_s']['obs'],
            DlnkWindow: self.sat_activity_params['min_duration_s']['dlnk'],
            XlnkWindow: self.sat_activity_params['min_duration_s']['xlnk']
        }

        # this is now less useful than I thought
        self.standard_activities = [ObsWindow,DlnkWindow,XlnkWindow]

        #  this should be as small as possible to prevent ill conditioning, but big enough that score factor constraints are still valid. 
        #  note: the size of this value is checked in make_model() below
        self.big_M_lat = 1e6


        # this big M should be conformatably bigger than the duration of any activity. It's used to switch on/off constraints that are related to differences in start/end times of activities, and in the worst case we expect such differences to max out at the legnth of the longest-duration activity. As long as we provide some comfortable margin past that duration, should be fine.
        # setting to 20 minutes
        self.big_M_act_t_dur_s = 30*60

        # this big M should be twice as large as any activity data volume 
        self.big_M_act_dv = 200000 # in Mb

        # see docs in super class
        self.allow_act_timing_constr_violations = False


    def run_sat_state_precheck(self,ecl_winds):
        """ propagate state forward from the initial state for every satellite to see if they're able to stay within resource limitations.  if we pass this check, then the MILP model should be solvable"""

        for sat_indx in range(self.num_sats):
            parsed_power_params = {
                "sat_edot_by_mode": self.sats_edot_by_mode_W[sat_indx],
                "sat_batt_storage": {'e_min': self.sats_emin_Wh[sat_indx],'e_max': self.sats_emax_Wh[sat_indx]},
                "power_units": None
            } 

            #  propagate energy storage from initial state through the whole planning window with no activities performed.  this is just checking if we can still meet any constraints while scheduling no activities whatsoever on this satellite.  in general, the overall planning system should not allow a satellite to get so low in energy usage that it might risk going below minimum without performing any activities whatsoever
            final_ES_no_acts,ES_state_went_below_min = propagate_sat_ES(
                self.planning_start_dt,
                self.planning_end_dt,
                sat_indx,
                self.sats_init_estate_Wh[sat_indx],
                [],
                ecl_winds[sat_indx],
                parsed_power_params,
                self.resource_delta_t_s
            )

            #  if the energy state did go below the minimum, we're not able to schedule with the satellite included (the MILP simply won't solve)
            # todo: what should we do if this actually happens? currently we're just assuming that it can't ... perhaps try scheduling with one of the satellites not included?
            if ES_state_went_below_min:
                raise RuntimeWarning('energy storage went below minimum for satellite %d in planning window (%s,%s) with initial energy storage state %f'%(sat_indx,self.planning_start_dt,self.planning_end_dt,self.sats_init_estate_Wh[sat_indx]))

    

    # def extract_resource_usage( self, decimation_factor =1, verbose = False):

    #     energy_usage = {}

    #     t_vals = []
    #     e_vals = [[] for sat_indx in range ( self.num_sats)]

    #     # note that this extraction uses the energy variables from the optimization, which are currently not constrained to be exactly equal to the energy delta from t-1 to t; they are merely bounded by it. Todo: extract concrete values based on activity execution times
    #     # TODO: this code feels super inefficient somehow.  make it better?
    #     for sat_indx, sat in enumerate (self.model.sats):
    #         last_tp_indx = 0
    #         for tp_indx in self.model.es_timepoint_indcs:

    #             if (tp_indx - last_tp_indx) < decimation_factor:
    #                 continue
    #             else:
    #                 e_vals[sat_indx].append(pe.value(self.model.var_sats_estore[sat,tp_indx]))

    #                 if sat_indx == 0:
    #                     t_vals.append(self.es_time_getter_dc.get_tp_from_tp_indx(tp_indx,out_units='minutes'))

    #     energy_usage['time_mins'] = t_vals
    #     energy_usage['e_sats'] = e_vals


    #     data_usage = {}

    #     t_vals = []
    #     d_vals = [[] for sat_indx in range ( self.num_sats)]

    #     # TODO: this code feels super inefficient somehow.  make it better?
    #     for sat_indx, sat in enumerate (self.model.sats):
    #         last_tp_indx = 0
    #         for tp_indx in self.model.ds_timepoint_indcs:

    #             if (tp_indx - last_tp_indx) < decimation_factor:
    #                 continue
    #             else:
    #                 d_vals[sat_indx].append(pe.value(self.model.var_sats_dstore[sat,tp_indx]))

    #                 if sat_indx == 0:
    #                     t_vals.append(self.ds_time_getter_dc.get_tp_from_tp_indx(tp_indx,out_units='minutes'))

    #     data_usage['time_mins'] = t_vals
    #     data_usage['d_sats'] = d_vals

    #     return  energy_usage, data_usage