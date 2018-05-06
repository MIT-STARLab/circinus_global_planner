# contains model and Solver for global planner activity scheduling capability
# 
# @author Kit Kennedy
#

from  datetime import timedelta
from copy import  deepcopy
from math import ceil

from pyomo import environ  as pe
from pyomo import opt  as po

import numpy as np

from circinus_tools  import time_tools as tt
from circinus_tools  import  constants as const
from circinus_tools  import io_tools
from circinus_tools.scheduling.custom_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow
from circinus_tools.scheduling.schedule_objects import Dancecard
from circinus_tools.scheduling.routing_objects import DataMultiRoute

class GPActivityScheduling():
    """docstring for GP activity scheduling"""
    
    # how close a binary variable must be to zero or one to be counted as that
    binary_epsilon = 0.1

    def __init__(self,gp_params):
        """initializes based on parameters
        
        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """

        scenario_params = gp_params['orbit_prop_params']['scenario_params']
        sat_params = gp_params['orbit_prop_params']['sat_params']
        as_params = gp_params['gp_general_params']['activity_scheduling_params']
        gp_inst_planning_params = gp_params['gp_instance_params']['planning_params']

        self.gp_agent_ID = gp_params['gp_instance_params']['gp_agent_ID']

        # records if the formulation model has been constructed successfully
        self.model_constructed = False

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
        # this is the mimimum latency requirement for the highest latency score factor, 1.0. If multiple routes/dlnks for a single obs have latency less than this, they will both have sf 1.0
        self.min_latency_for_sf_1_mins =as_params['min_latency_for_sf_1_mins']

        self.planning_start_dt  = gp_inst_planning_params['planning_start_dt']
        self.planning_end_obs_xlnk_dt = gp_inst_planning_params['planning_end_obs_xlnk_dt']
        self.planning_end_dlnk_dt  = gp_inst_planning_params['planning_end_dlnk_dt']
        self.planning_end_dt  = self.planning_end_dlnk_dt

        #  the "effectively zero" number.
        self.dv_epsilon = as_params['dv_epsilon_Mb']
        self.resource_margin_obj_num_timepoints = as_params['resource_margin_obj_num_timepoints']

        self.obj_weights =as_params['obj_weights']

        self.use_symmetric_xlnk_windows = gp_params['gp_general_params']['other_params']['use_symmetric_xlnk_windows']

        self.power_params = sat_params['power_params_sorted']
        self.data_storage_params = sat_params['data_storage_params_sorted']
        self.initial_state = sat_params['initial_state_sorted']
        sat_activity_params = sat_params['activity_params']


        self.sats_dmin_Mb = [1000*ds_params['data_storage_Gbit']['d_min'][ds_params['storage_option']] for ds_params in self.data_storage_params]
        self.sats_dmax_Mb = [1000*ds_params['data_storage_Gbit']['d_max'][ds_params['storage_option']] for ds_params in self.data_storage_params]

        self.energy_unit = "Wh"

        # these lists are in order of satellite index because we've sorted 
        self.sats_init_estate_Wh = [sat_state['batt_e_Wh'] for sat_state in self.initial_state]
        self.sats_edot_by_mode_W = []
        self.sats_emin_Wh = []
        self.sats_emax_Wh = []
        for p_params in self.power_params:
            sat_edot_by_mode,sat_batt_storage = io_tools.parse_power_consumption_params(p_params)

            self.sats_edot_by_mode_W.append (sat_edot_by_mode)
            self.sats_emin_Wh.append (sat_batt_storage['e_min'])
            self.sats_emax_Wh.append (sat_batt_storage['e_max'])

        self.act_transition_time_map = {
            ('inter-sat',DlnkWindow,DlnkWindow): sat_activity_params['transition_time_s']['inter-sat']['dlnk-dlnk'],
            "default":  sat_activity_params['transition_time_s']['default']
        }

        self.min_act_duration_s = {
            ObsWindow: sat_activity_params['min_duration_s']['obs'],
            DlnkWindow: sat_activity_params['min_duration_s']['dlnk'],
            XlnkWindow: sat_activity_params['min_duration_s']['xlnk']
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

    def gen_inter_act_constraint(self,var_list,constr_list,transition_time_req,model_objs_act1,model_objs_act2):

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

                # If they don't overlap in center time, but they do overlap to some amount, then we need to constrain their end and start times to be consistent with one another
                else:
                    M = max(act1.duration,act2.duration).total_seconds()
                    constr_disable_1 = M*(1-model_objs_act1['var_act_indic'])
                    constr_disable_2 = M*(1-model_objs_act2['var_act_indic'])
        
                    constr = center_time_diff - time_adjust_1 - time_adjust_2 + constr_disable_1 + constr_disable_2 >= transition_time_req

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

            return constr, var_constr_violation, min_constr_violation


    def gen_intra_sat_act_overlap_constraints(self,c_overlap,c_duration,sats_acts,act_model_objs_getter,constraint_violation_model_objs):

        # keep track of this, so we know to warn user about default transition time usage
        #  note: this is not a parameter! it's merely helpful for reporting warnings
        used_default_transition_time = False

        intra_sat_act_constr_violation_acts_list = constraint_violation_model_objs['intra_sat_act_constr_violation_acts_list']
        var_intra_sat_act_constr_violations = constraint_violation_model_objs['var_intra_sat_act_constr_violations']
        intra_sat_act_constr_bounds = constraint_violation_model_objs['intra_sat_act_constr_bounds']
        min_var_intra_sat_act_constr_violation_list = constraint_violation_model_objs['min_var_intra_sat_act_constr_violation_list']

        for sat_indx in range (self.num_sats):
            num_sat_acts = len(sats_acts[sat_indx])
            for  first_act_indx in  range (num_sat_acts):
                act1 = sats_acts[sat_indx][first_act_indx]
                # note: generally a function from the AS subclass
                model_objs_act1 = act_model_objs_getter(act1)

                length_1 = model_objs_act1['var_dv_utilization']/act1.ave_data_rate
                c_duration.add( length_1 >= model_objs_act1['var_act_indic'] * self.min_act_duration_s[type(act1)])
                
                for  second_act_indx in  range (first_act_indx+1,num_sat_acts):
                    act2 = sats_acts[sat_indx][second_act_indx]

                    # act list should be sorted
                    assert(act2.center >= act1.center)

                    # get the transition time requirement between these activities
                    try:
                        transition_time_req = self.act_transition_time_map[("intra-sat",type(act1),type(act2))]
                    # if not explicitly specified, go with default transition time requirement
                    except KeyError:
                        used_default_transition_time = True
                        transition_time_req = self.act_transition_time_map["default"]

                    # if there is enough transition time between the two activities, no constraint needs to be added
                    #  note that we are okay even if for some reason Act 2 starts before Act 1 ends, because time deltas return negative total seconds as well
                    if (act2.start - act1.end).total_seconds() >= transition_time_req:
                        #  don't need to do anything,  continue on to next activity pair
                        continue

                    else:
                        # note: generally a function from the AS subclass
                        model_objs_act2 = act_model_objs_getter(act2)

                        constr, var_constr_violation, min_constr_violation = self.gen_inter_act_constraint(
                            var_intra_sat_act_constr_violations,
                            intra_sat_act_constr_bounds,
                            transition_time_req,
                            model_objs_act1,
                            model_objs_act2
                        )

                        #  add the constraint, regardless of whether or not it's a "big M" constraint, or a constraint violation constraint - they're handled the same
                        c_overlap.add( constr )

                        #  if it's a constraint violation constraint, then we have a variable to deal with
                        if not min_constr_violation is None:
                            min_var_intra_sat_act_constr_violation_list.append(min_constr_violation)
                            intra_sat_act_constr_violation_acts_list.append((act1,act2))

        return used_default_transition_time


    def gen_inter_sat_act_overlap_constraints(self,c_overlap,sats_dlnks,act_model_objs_getter,constraint_violation_model_objs):

        # keep track of this, so we know to warn user about default transition time usage
        #  note: this is not a parameter! it's merely helpful for reporting warnings
        used_default_transition_time = False

        inter_sat_act_constr_violation_acts_list = constraint_violation_model_objs['inter_sat_act_constr_violation_acts_list']
        var_inter_sat_act_constr_violations = constraint_violation_model_objs['var_inter_sat_act_constr_violations']
        inter_sat_act_constr_bounds = constraint_violation_model_objs['inter_sat_act_constr_bounds']
        min_var_inter_sat_act_constr_violation_list = constraint_violation_model_objs['min_var_inter_sat_act_constr_violation_list']

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
                        try:
                            transition_time_req = self.act_transition_time_map[("inter-sat",DlnkWindow,DlnkWindow)]
                        # if not explicitly specified, go with default transition time requirement
                        except KeyError:
                            used_default_transition_time = True
                            transition_time_req = self.act_transition_time_map["default"]

                        # if there is enough transition time between the two activities, no constraint needs to be added
                        #  note that we are okay even if for some reason Act 2 starts before Act 1 ends, because time deltas return negative total seconds as well
                        if (act2.start - act1.end).total_seconds() >= transition_time_req:
                            #  don't need to do anything,  continue on to next activity pair
                            continue

                        else:
                            # note: generally a function from the AS subclass
                            model_objs_act1 = act_model_objs_getter(act1)
                            # note: generally a function from the AS subclass
                            model_objs_act2 = act_model_objs_getter(act2)
                        
                            constr, var_constr_violation, min_constr_violation = self.gen_inter_act_constraint(
                                var_inter_sat_act_constr_violations,
                                inter_sat_act_constr_bounds,
                                transition_time_req,
                                model_objs_act1,
                                model_objs_act2
                            )

                            #  add the constraint, regardless of whether or not it's a "big M" constraint, or a constraint violation constraint - they're handled the same
                            c_overlap.add( constr )

                            #  if it's a constraint violation constraint, then we have a variable to deal with
                            if not min_constr_violation is None:
                                # model.var_inter_sat_act_constr_violations.add(var_constr_violation)
                                min_var_inter_sat_act_constr_violation_list.append(min_constr_violation)
                                inter_sat_act_constr_violation_acts_list.append((act1,act2))

        return used_default_transition_time

    def solve(self):

        if not self.model_constructed:
            return

        solver = po.SolverFactory(self.solver_name)
        if self.solver_name == 'gurobi':
            # note default for this is 1e-4, or 0.01%
            solver.options['TimeLimit'] = self.solver_params['max_runtime_s']
            solver.options['MIPGap'] = self.solver_params['optimality_gap']
            solver.options['IntFeasTol'] = self.solver_params['integer_feasibility_tolerance']
            # other options...
            # solver.options['Cuts'] = 0
            # solver.options['MIPFocus'] = 1 #for finding feasible solutions quickly
            # solver.options['MIPFocus'] = 3 #for lowering the mip gap

        elif self.solver_name == 'cplex':
            solver.options['timelimit'] = self.solver_params['max_runtime_s']
            solver.options['mip_tolerances_mipgap'] = self.solver_params['optimality_gap']
            solver.options['mip_tolerances_integrality'] = self.solver_params['integer_feasibility_tolerance']


        # if we're running things remotely, then we will use the NEOS server (https://neos-server.org/neos/)
        if self.solver_params['run_remotely']:
            solver_manager = po.SolverManagerFactory('neos')
            results = solver_manager.solve(self.model, opt= solver)
        else:
            # tee=True displays solver output in the terminal
            # keepfiles=True  keeps files passed to and from the solver
            results =  solver.solve(self.model, tee=True, keepfiles= False)

        if (results.solver.status == po.SolverStatus.ok) and (results.solver.termination_condition == po.TerminationCondition.optimal):
            print('this is feasible and optimal')
        elif results.solver.termination_condition == po.TerminationCondition.infeasible:
            print ('infeasible')
            raise RuntimeError('Model is infeasible with current parameters')
        else:
            # something else is wrong
            print (results.solver)

    def print_sol_all(self):
        for v in self.model.component_objects(pe.Var, active=True):
            print ("Variable",v)
            varobject = getattr(self.model, str(v))
            for index in varobject:
                print (" ",index, varobject[index].value)

    def extract_resource_usage( self, decimation_factor =1, verbose = False):

        energy_usage = {}

        t_vals = []
        e_vals = [[] for sat_indx in range ( self.num_sats)]

        # note that this extraction uses the energy variables from the optimization, which are currently not constrained to be exactly equal to the energy delta from t-1 to t; they are merely bounded by it. Todo: extract concrete values based on activity execution times
        # TODO: this code feels super inefficient somehow.  make it better?
        for sat_indx, sat in enumerate (self.model.sats):
            last_tp_indx = 0
            for tp_indx in self.model.es_timepoint_indcs:

                if (tp_indx - last_tp_indx) < decimation_factor:
                    continue
                else:
                    e_vals[sat_indx].append(pe.value(self.model.var_sats_estore[sat,tp_indx]))

                    if sat_indx == 0:
                        t_vals.append(self.es_time_getter_dc.get_tp_from_tp_indx(tp_indx,out_units='minutes'))

        energy_usage['time_mins'] = t_vals
        energy_usage['e_sats'] = e_vals


        data_usage = {}

        t_vals = []
        d_vals = [[] for sat_indx in range ( self.num_sats)]

        # TODO: this code feels super inefficient somehow.  make it better?
        for sat_indx, sat in enumerate (self.model.sats):
            last_tp_indx = 0
            for tp_indx in self.model.ds_timepoint_indcs:

                if (tp_indx - last_tp_indx) < decimation_factor:
                    continue
                else:
                    d_vals[sat_indx].append(pe.value(self.model.var_sats_dstore[sat,tp_indx]))

                    if sat_indx == 0:
                        t_vals.append(self.ds_time_getter_dc.get_tp_from_tp_indx(tp_indx,out_units='minutes'))

        data_usage['time_mins'] = t_vals
        data_usage['d_sats'] = d_vals

        return  energy_usage, data_usage