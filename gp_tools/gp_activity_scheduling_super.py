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
from circinus_tools.scheduling.formulation.schedulers import PyomoMILPScheduling

class GPActivityScheduling(PyomoMILPScheduling):
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
        
        self.latency_params = gp_params['gp_general_params']['other_params']['latency_calculation']
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
        self.epsilon_fixed_utilization = as_params['epsilon_fixed_utilization']
        # this is the mimimum latency requirement for the highest latency score factor, 1.0. If multiple routes/dlnks for a single obs have latency less than this, they will both have sf 1.0
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

        # allow activities to overlap, and penalize them for doing so. The code should work, but hasn't been extensively vetted for its usefulness (does seem surprisingly unresponsive to changing weights for constraint violation in the obj function...). Note having these violations allowed generally won't play well with extracting routes in coupled AS due to data route validation checks. 
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

    def gen_inter_act_constraint(self,var_list,constr_list,transition_time_req,model_objs_act1,model_objs_act2):
            #  the regular constraint is the constraint that is enforced in the mixed integer linear program
            #  the "binding expression"  is the expression that, when evaluated after MILP solution, will be zero if the constraint is binding and greater than zero if the constraint is not binding ( i.e., if the constraint is not binding, there is room for one of the variables to grow). The constraint should have any big M factors removed in the binding expression, because those obfuscate information about whether the constraint is binding or not. The binding expression shall may be negative if big M components are normally present.

            var_constr_violation = None
            min_constr_violation = None

            act1 = model_objs_act1['act_object']
            act2 = model_objs_act2['act_object']

            center_time_diff = (act2.center - act1.center).total_seconds()
            # var is a variable, the amount of dv used for the link
            time_adjust_1 = model_objs_act1['var_dv_utilization']/2/act1.ave_data_rate
            time_adjust_2 = model_objs_act2['var_dv_utilization']/2/act2.ave_data_rate
            # par is a parameter
            max_time_adjust_1 = model_objs_act1['par_dv_capacity']*1/2/act1.ave_data_rate
            max_time_adjust_2 = model_objs_act2['par_dv_capacity']*1/2/act2.ave_data_rate

            if not self.allow_act_timing_constr_violations:
                #  if the activities overlap in center time (including  transition time), then it's not possible to have sufficient transition time between them.  only allow one
                if (act2.center - act1.center).total_seconds() <= transition_time_req:
                    constr = model_objs_act1['var_act_indic']+ model_objs_act2['var_act_indic'] <= 1
                    binding_expr = 1 - model_objs_act1['var_act_indic'] - model_objs_act2['var_act_indic']

                # If they don't overlap in center time, but they do overlap to some amount, then we need to constrain their end and start times to be consistent with one another
                else:
                    M = max(act1.duration,act2.duration).total_seconds()
                    constr_disable_1 = M*(1-model_objs_act1['var_act_indic'])
                    constr_disable_2 = M*(1-model_objs_act2['var_act_indic'])
        
                    constr = center_time_diff - time_adjust_1 - time_adjust_2 + constr_disable_1 + constr_disable_2 >= transition_time_req
                    binding_expr = center_time_diff - time_adjust_1 - time_adjust_2 - transition_time_req 

            else:
                # time_adjust_N can go as low as zero, so the constraint violation can be this at its lowest
                min_constr_violation = center_time_diff - max_time_adjust_1 - max_time_adjust_2 - transition_time_req

                assert(min_constr_violation < 0)

                # deal with adding a variable to represent this constraint violation. ( I hate that I have to do it this way, deep within this function, but it seems like the approach you have to use for dynamically generating variable lists in pyomo, bizarrely. refer to: https://projects.coin-or.org/Coopr/browser/pyomo/trunk/pyomo/core/base/var.py?rev=11067, https://groups.google.com/forum/#!topic/pyomo-forum/modS1VkPxW0
                var_list.add()
                var_constr_violation = var_list[len(var_list)]

                #  bounds on range of constraint violation variable

                # bound minimum of the constraint violation (where both activities have 1.0 utilization), to keep the problem formulation as tight as possible
                var_constr_violation.setlb(min_constr_violation)
                # want the constraint violation only to go to zero at its maximum, because we don't want to reward times where there is no constraint violation, only penalize
                var_constr_violation.setub(0)

                #  the actual time constraint that bounds the constraint violation
                constr = center_time_diff - time_adjust_1 - time_adjust_2 - transition_time_req >= var_constr_violation
                binding_expr = center_time_diff - time_adjust_1 - time_adjust_2 - transition_time_req - var_constr_violation

            return constr, binding_expr,var_constr_violation, min_constr_violation


    def gen_intra_sat_act_overlap_constraints(self,c_overlap,c_duration,sats_acts,act_model_objs_getter,constraint_violation_model_objs):

        intra_sat_act_constr_violation_acts_list = constraint_violation_model_objs['intra_sat_act_constr_violation_acts_list']
        var_intra_sat_act_constr_violations = constraint_violation_model_objs['var_intra_sat_act_constr_violations']
        intra_sat_act_constr_bounds = constraint_violation_model_objs['intra_sat_act_constr_bounds']
        min_var_intra_sat_act_constr_violation_list = constraint_violation_model_objs['min_var_intra_sat_act_constr_violation_list']

        binding_expr_overlap_by_act = {}
        binding_expr_duration_by_act = {}

        for sat_indx in range (self.num_sats):
            num_sat_acts = len(sats_acts[sat_indx])
            for  first_act_indx in  range (num_sat_acts):
                act1 = sats_acts[sat_indx][first_act_indx]
                # note: generally a function from the AS subclass
                model_objs_act1 = act_model_objs_getter(act1)

                length_1 = model_objs_act1['var_dv_utilization']/act1.ave_data_rate
                c_duration.add( length_1 >= model_objs_act1['var_act_indic'] * self.min_act_duration_s[type(act1)])
                binding_expr_duration_by_act[act1] = length_1 - model_objs_act1['var_act_indic'] * self.min_act_duration_s[type(act1)]
                
                for  second_act_indx in  range (first_act_indx+1,num_sat_acts):
                    act2 = sats_acts[sat_indx][second_act_indx]

                    # act list should be sorted
                    assert(act2.center >= act1.center)

                    # get the transition time requirement between these activities
                    transition_time_req = io_tools.get_transition_time_req(act1,act2,sat_indx,sat_indx,self.sat_activity_params)

                    # if there is enough transition time between the two activities, no constraint needs to be added
                    #  note that we are okay even if for some reason Act 2 starts before Act 1 ends, because time deltas return negative total seconds as well
                    if (act2.original_start - act1.original_end).total_seconds() >= transition_time_req:
                        #  don't need to do anything,  continue on to next activity pair
                        continue

                    else:
                        # note: generally a function from the AS subclass
                        model_objs_act2 = act_model_objs_getter(act2)

                        constr, binding_expr, var_constr_violation, min_constr_violation = self.gen_inter_act_constraint(
                            var_intra_sat_act_constr_violations,
                            intra_sat_act_constr_bounds,
                            transition_time_req,
                            model_objs_act1,
                            model_objs_act2
                        )

                        #  add the constraint, regardless of whether or not it's a "big M" constraint, or a constraint violation constraint - they're handled the same
                        c_overlap.add( constr )
                        binding_expr_overlap_by_act.setdefault(act1,[]).append(binding_expr)
                        binding_expr_overlap_by_act.setdefault(act2,[]).append(binding_expr)

                        #  if it's a constraint violation constraint, then we have a variable to deal with
                        if not min_constr_violation is None:
                            min_var_intra_sat_act_constr_violation_list.append(min_constr_violation)
                            intra_sat_act_constr_violation_acts_list.append((act1,act2))

        return binding_expr_duration_by_act,binding_expr_overlap_by_act


    def gen_inter_sat_act_overlap_constraints(self,c_overlap,sats_dlnks,act_model_objs_getter,constraint_violation_model_objs):

        inter_sat_act_constr_violation_acts_list = constraint_violation_model_objs['inter_sat_act_constr_violation_acts_list']
        var_inter_sat_act_constr_violations = constraint_violation_model_objs['var_inter_sat_act_constr_violations']
        inter_sat_act_constr_bounds = constraint_violation_model_objs['inter_sat_act_constr_bounds']
        min_var_inter_sat_act_constr_violation_list = constraint_violation_model_objs['min_var_inter_sat_act_constr_violation_list']

        binding_expr_overlap_by_act = {}

        for sat_indx in range (self.num_sats):
            num_sat_acts = len(sats_dlnks[sat_indx])
            
            for other_sat_indx in range (self.num_sats):
                if other_sat_indx == sat_indx:
                    continue

                num_other_sat_acts = len(sats_dlnks[other_sat_indx])

                for  sat_act_indx in  range (num_sat_acts):

                    act1 = sats_dlnks[sat_indx][sat_act_indx]
                    
                    for  other_sat_act_indx in  range (num_other_sat_acts):
                        act2 = sats_dlnks[other_sat_indx][other_sat_act_indx]

                        assert(type(act1) == DlnkWindow and type(act2) == DlnkWindow)

                        # this line is pretty important - only consider overlap if they're looking at the same GS. I forgot to add this before and spent days wondering why the optimization process was progressing so slowly (hint: it's really freaking constrained and there's not much guidance for finding a good objective value if no downlink can overlap in time with any other downlink)
                        if act1.gs_indx != act2.gs_indx:
                            continue

                        # we're considering windows across satellites, so they could be out of order temporally. These constraints are only valid if act2 is after act1 (center time). Don't worry, as we loop through satellites, we consider both directions (i.e. act1 and act2 will be swapped in another iteration, and we'll get past this check and impose the required constraints)
                        if (act2.center - act1.center).total_seconds() < 0:
                            continue

                        # get the transition time requirement between these activities
                        transition_time_req = io_tools.get_transition_time_req(act1,act2,sat_indx,other_sat_indx,self.sat_activity_params)                   

                        # if there is enough transition time between the two activities, no constraint needs to be added
                        #  note that we are okay even if for some reason Act 2 starts before Act 1 ends, because time deltas return negative total seconds as well
                        if (act2.original_start - act1.original_end).total_seconds() >= transition_time_req:
                            #  don't need to do anything,  continue on to next activity pair
                            continue

                        else:
                            # note: generally a function from the AS subclass
                            model_objs_act1 = act_model_objs_getter(act1)
                            # note: generally a function from the AS subclass
                            model_objs_act2 = act_model_objs_getter(act2)
                        
                            constr, binding_expr, var_constr_violation, min_constr_violation = self.gen_inter_act_constraint(
                                var_inter_sat_act_constr_violations,
                                inter_sat_act_constr_bounds,
                                transition_time_req,
                                model_objs_act1,
                                model_objs_act2
                            )

                            #  add the constraint, regardless of whether or not it's a "big M" constraint, or a constraint violation constraint - they're handled the same
                            c_overlap.add( constr )
                            binding_expr_overlap_by_act.setdefault(act1,[]).append(binding_expr)
                            binding_expr_overlap_by_act.setdefault(act2,[]).append(binding_expr)

                            #  if it's a constraint violation constraint, then we have a variable to deal with
                            if not min_constr_violation is None:
                                # model.var_inter_sat_act_constr_violations.add(var_constr_violation)
                                min_var_inter_sat_act_constr_violation_list.append(min_constr_violation)
                                inter_sat_act_constr_violation_acts_list.append((act1,act2))

        return binding_expr_overlap_by_act

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