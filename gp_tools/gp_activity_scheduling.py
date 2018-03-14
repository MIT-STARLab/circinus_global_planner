# Contains functionality for turning input data structures into the 
# objects used by the global planner  scheduling module
# 
# 
# @author Kit Kennedy
#
#  note that a path is the same as a route. 

from  datetime import timedelta

from pyomo import environ  as pe
from pyomo import opt  as po

from .custom_activity_window import   ObsWindow,  DlnkWindow, XlnkWindow
from .routing_objects import DataRoute
from .schedule_objects import Dancecard
from circinus_tools  import io_tools

class GPActivityScheduling():
    """docstring for GP activity scheduling"""
    
    # how close a binary variable must be to zero or one to be counted as that
    binary_epsilon = 0.1

    def __init__(self,params):
        self.solver_max_runtime =params['solver_max_runtime_s']
        self.solver_name =params['solver_name']
        self.solver_run_remotely =params['solver_run_remotely']
        self.min_path_dv =params['min_path_dv_Mb']
        self.num_sats=params['num_sats']
        self.transition_time_s=params['transition_time_s']
        self.start_utc_dt  =params['start_utc_dt']
        self.end_utc_dt  =params['end_utc_dt']
        self.resource_delta_t_s  =params['resource_delta_t_s']
        sat_id_order= params['sat_id_order']

        #  sort these now in case they weren't before
        self.power_params = io_tools.sort_input_params_by_sat_indcs(params['power_params'],sat_id_order)
        self.initial_state = io_tools.sort_input_params_by_sat_indcs(params['initial_state'],sat_id_order)

        # these lists are in order of satellite index because we've sorted 
        self.sats_init_estate_Wh = [sat_state['batt_e_Wh'] for sat_state in self.initial_state]
        self.sats_emin_Wh = [p_params['battery_storage_Wh']['e_min'][p_params['battery_option']] for p_params in self.power_params]
        self.sats_emax_Wh = [p_params['battery_storage_Wh']['e_max'][p_params['battery_option']] for p_params in self.power_params]

        self.sats_edot_by_act_W = []
        for p_params in self.power_params:
            sat_edot_by_act = {}
            sat_edot_by_act['base'] = p_params['power_consumption_W']['base'][p_params['base_option']]
            sat_edot_by_act['obs'] = p_params['power_consumption_W']['obs'][p_params['obs_option']]
            sat_edot_by_act['dlnk'] = p_params['power_consumption_W']['dlnk'][p_params['dlnk_option']]
            sat_edot_by_act['xlnk'] = p_params['power_consumption_W']['xlnk'][p_params['xlnk_option']]
            sat_edot_by_act['charging'] = p_params['power_consumption_W']['orbit_insunlight_average_charging'][p_params['charging_option']]
            self.sats_edot_by_act_W.append (sat_edot_by_act)

        self.act_type_map = {
            ObsWindow: 'obs',
            DlnkWindow: 'dlnk',
            XlnkWindow: 'xlnk'
        }

    def get_stats(self,verbose=True):
        stats = {}
        stats['num_acts'] = sum([len ( self.all_acts_indcs)])
        stats['num_obs_acts'] = sum([len ( self.obs_act_indcs)])
        stats['num_link_acts'] = sum([len ( self.link_act_indcs)])
        stats['num_variables'] = self.model.nvariables ()
        stats['num_nobjectives'] = self.model.nobjectives ()
        stats['num_nconstraints'] = self.model.nconstraints ()

        if verbose:
            print ( "Considering %d activity windows" % (stats['num_acts']))
            print ( "Considering %d observation windows" % (stats['num_obs_acts']))
            print ( "Considering %d link windows" % (stats['num_link_acts']))
            print ( 'self.model.nvariables ()')
            print ( self.model.nvariables ())
            print ( 'self.model.nobjectives ()')
            print ( self.model.nobjectives ())
            print ( 'self.model.nconstraints ()')
            print ( self.model.nconstraints ())

        return stats

    def get_activity_structs( self,routes_flat):

        #  all activities are uniquely indexed. these structures keep track of those, and the mapping to activity objects
        all_acts_indcs = []
        #  these structures are for lookup in both directions
        all_acts_by_indx = {}
        all_acts_by_obj = {}

        # these structures keep track of the subset of unique indices that correspond to observations and links. Also we keep track of what data routes correspond to an activity
        link_act_indcs = []
        obs_act_indcs = []
        path_indcs_by_link_act = {}
        path_indcs_by_obs_act = {}
        path_indcs_by_act = {}
        dv_by_link_act = {}
        dv_by_obs_act = {}
        dv_by_act = {}

        sat_acts = [[] for sat_indx in range (self.num_sats)]

        new_act_indx = 0
        for dr_indx, dr in enumerate (routes_flat):
            for act in dr.route:

                # if we haven't yet seen this activity, then add it to bookkeeping
                if not act in all_acts_by_obj.keys():
                    act_indx = new_act_indx
                    all_acts_indcs.append(act_indx)
                    all_acts_by_indx[act_indx] = act
                    all_acts_by_obj[act] = act_indx
                    path_indcs_by_act[act_indx] = []
                    path_indcs_by_act[act_indx].append (dr_indx)
                    dv_by_act[act_indx] = act.data_vol
                    new_act_indx += 1

                    # also need to add it to the list and dictionary for observations
                    if type(act) == ObsWindow:
                        sat_acts[act.sat_indx].append(act)
                        obs_act_indcs.append(act_indx)
                        path_indcs_by_obs_act[act_indx] = []
                        path_indcs_by_obs_act[act_indx].append (dr_indx)
                        dv_by_obs_act[act_indx] = act.data_vol

                    # also need to add it to the list and dictionary for links
                    if type(act) == DlnkWindow:
                        link_act_indcs.append(act_indx)
                        path_indcs_by_link_act[act_indx] = []
                        path_indcs_by_link_act[act_indx].append (dr_indx)
                        dv_by_link_act[act_indx] = act.data_vol
                        sat_acts[act.sat_indx].append(act)

                    if type(act) == XlnkWindow:
                        link_act_indcs.append(act_indx)
                        path_indcs_by_link_act[act_indx] = []
                        path_indcs_by_link_act[act_indx].append (dr_indx)
                        dv_by_link_act[act_indx] = act.data_vol
                        sat_acts[act.sat_indx].append(act)
                        sat_acts[act.xsat_indx].append(act)

                #  if we have already seen the activity,  then just need to update the appropriate structures
                else:
                    act_indx = all_acts_by_obj[act]
                    path_indcs_by_act[act_indx].append (dr_indx)

                    # add the data route index
                    if type(act) == ObsWindow:
                        path_indcs_by_obs_act[act_indx].append (dr_indx)
                    if type(act) == DlnkWindow or type(act) == XlnkWindow:
                        path_indcs_by_link_act[act_indx].append (dr_indx)

        # print (path_indcs_by_link_act)
        # print (path_indcs_by_obs_act)

        #  sort the activities, because we'll need that for constructing constraints
        for sat_indx in range (self.num_sats):
            sat_acts[sat_indx].sort(key=lambda x: x.start)

        return sat_acts,all_acts_indcs,path_indcs_by_act,dv_by_act,all_acts_by_obj,obs_act_indcs,path_indcs_by_obs_act,dv_by_obs_act,link_act_indcs,path_indcs_by_link_act,dv_by_link_act
                    




    def make_model ( self,routes_flat, verbose = True):
        model = pe.ConcreteModel()

        # print(routes_flat)

        self.routes_flat = routes_flat

        if verbose:
            pass

        ##############################
        #  Make indices/ subscripts
        ##############################

        # note: a path is the same as a route

        (sat_acts,
            all_acts_indcs,
            path_indcs_by_act,
            dv_by_act,
            all_acts_by_obj,
            obs_act_indcs,
            path_indcs_by_obs_act,
            dv_by_obs_act,
            link_act_indcs,
            path_indcs_by_link_act,
            dv_by_link_act) =  self.get_activity_structs(routes_flat)

        #  calculate center times and average data rates for activities in advance
        for sat_indx in range (self.num_sats):
            for act in sat_acts[sat_indx]:
                act.center_time = ( act.end - act.start) / 2

                act.ave_data_rate =  act.data_vol / ( act.end - act.start).total_seconds ()

        # construct a set of dance cards for every satellite, 
        # each of which keeps track of all of the activities of satellite 
        # can possibly execute at any given time slice delta T. 
        activity_dancecards = [Dancecard(self.start_utc_dt,self.end_utc_dt,self.resource_delta_t_s) for sat_indx in range (self.num_sats)]
        for sat_indx in range (self.num_sats): 
            activity_dancecards[sat_indx].add_winds_to_dancecard(sat_acts[sat_indx])

        self.all_acts_indcs = all_acts_indcs
        self.obs_act_indcs = obs_act_indcs
        self.link_act_indcs = link_act_indcs

        #  subscript for each path p
        model.paths = pe.Set(initialize= range(len(routes_flat)))
        #  subscript for each activity a
        model.acts = pe.Set(initialize= all_acts_indcs)
        #  subscript for each satellite
        model.sats = pe.Set(initialize=  range ( self.num_sats))

        # timepoints is the indices, whereas timepoints_s is the time values in seconds
        #  NOTE: we assume the same time system for every satellite
        model.timepoints = pe.Set(initialize=  activity_dancecards[0].get_timepoint_indices ())
        self.timepoints_s = activity_dancecards[0].get_timepoint_values(units='seconds')
        self.timepoints_m = activity_dancecards[0].get_timepoint_values(units='minutes')

        #  unique indices for observation and link acts
        model.obs_acts = pe.Set(initialize= obs_act_indcs)
        model.link_acts = pe.Set(initialize= link_act_indcs)

        ##############################
        #  Make parameters
        ##############################

        model.par_min_path_dv = pe.Param (initialize=self.min_path_dv)
        model.par_obs_dv = pe.Param(model.obs_acts,initialize =dv_by_obs_act)
        model.par_link_dv = pe.Param(model.link_acts,initialize =dv_by_link_act)
        model.par_act_dv = pe.Param(model.acts,initialize =dv_by_act)
        model.par_path_dv = pe.Param(model.paths,initialize ={ i: dr.data_vol for i,dr in enumerate (routes_flat)})
        # each of these is essentially a dictionary indexed by link or obs act indx, with  the value being a list of path indices that are included within that act
        # these are valid indices into model.paths
        model.par_path_subscrs_by_link_act = pe.Param(model.link_acts,initialize =path_indcs_by_link_act)
        model.par_path_subscrs_by_obs_act = pe.Param(model.obs_acts,initialize =path_indcs_by_obs_act)
        model.par_path_subscrs_by_act = pe.Param(model.acts,initialize =path_indcs_by_act)

        model.par_resource_delta_t = pe.Param (initialize= self.resource_delta_t_s) 
        model.par_sats_estore_initial = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_init_estate_Wh)})
        model.par_sats_estore_min = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_emin_Wh)})
        model.par_sats_estore_max = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_emax_Wh)})
        model.par_sats_edot_by_act = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_edot_by_act_W)})


        ##############################
        #  Make variables
        ##############################


        # activity utilization variable indicating how much of an activity's capacity is used [1]
        model.var_activity_utilization  = pe.Var (model.acts, bounds =(0,1))
        # path utilization variable indicating how much of a path's capacity is used [2]
        model.var_path_utilization  = pe.Var (model.paths, bounds =(0,1))
        #  indicator variables for whether or not paths [3] and activities [4] have been chosen
        model.var_path_indic  = pe.Var (model.paths, within = pe.Binary)
        model.var_act_indic  = pe.Var (model.acts, within = pe.Binary)
        
        # satellite energy storage
        model.var_sats_estore  = pe.Var (model.sats,  model.timepoints,  within = pe.NonNegativeReals)
        

        ##############################
        #  Make constraints
        ##############################

        # TODO: renumber  these with the final numbering

        def c1_rule( model,a):
            return (model.par_act_dv[a]*model.var_activity_utilization[a] -
                        sum(model.par_path_dv[p]*model.var_path_utilization[p] 
                            for p in model.par_path_subscrs_by_act[a]) 
                    >= 0)
        model.c1 =pe.Constraint ( model.acts,  rule=c1_rule)

        def c2_rule( model,p):
            return model.par_path_dv[p]*model.var_path_utilization[p] >= model.par_min_path_dv*model.var_path_indic[p]
        model.c2 =pe.Constraint ( model.paths,  rule=c2_rule)

        def c3_rule( model,a):
            return model.var_act_indic[a] >=  model.var_activity_utilization[a]
        model.c3 =pe.Constraint ( model.acts,  rule=c3_rule)        

        #  activity overlap constraints [4],[5]
        model.c4  = pe.ConstraintList()
        model.c5  = pe.ConstraintList()
        for sat_indx in range (self.num_sats):
            num_sat_acts = len(sat_acts[sat_indx])
            for  first_act_indx in  range (num_sat_acts):
                for  second_act_indx in  range (first_act_indx+1,num_sat_acts):
                    act1 = sat_acts[sat_indx][first_act_indx]
                    act2 = sat_acts[sat_indx][ second_act_indx]
                    # get the unique indices into model.acts
                    act1_uindx = all_acts_by_obj[act1]
                    act2_uindx = all_acts_by_obj[act2]

                    # if there is enough transition time between the two activities, no constraint needs to be added
                    if (act2.start - act1.end).total_seconds() >= self.transition_time_s['simple']:
                        #  don't need to do anything,  continue on to next activity pair
                        continue

                    #  if the activities overlap in center time, then it's not possible to have sufficient transition time between them
                    #  add constraint to rule out the possibility of scheduling both of them
                    elif (act2.center_time - act1.center_time).total_seconds() <= self.transition_time_s['simple']:
                        model.c4.add( model.var_act_indic[act1_uindx]+ model.var_act_indic[act2_uindx] <= 1)

                    # If they don't overlap in center time, but they do overlap to some amount, then we need to constrain their end and start times to be consistent with one another
                    else:
                        center_time_diff = (act2.center_time - act1.center_time).total_seconds()
                        # this is the adjustment added to the center time to get to the start or end of the activity
                        time_adjust_1 = model.par_act_dv[act1_uindx]*model.var_activity_utilization[act1_uindx]/2/act1.ave_data_rate
                        time_adjust_2 = model.par_act_dv[act2_uindx]*model.var_activity_utilization[act2_uindx]/2/act2.ave_data_rate
                        model.c5.add( center_time_diff - time_adjust_1 - time_adjust_2 >= self.transition_time_s['simple'])

        #  energy constraints [6]
        model.c6  = pe.ConstraintList()
        for sat_indx in range (self.num_sats): 

            # timepoint serves as an index into the satellite activity dance cards
            # we will use timepoint as an index into model data structures for consistency
            for tp_indx, timepoint in enumerate(model.timepoints):
                #  constraining first time step to initial energy storage
                #  continue for loop afterwards because no geq/leq constraints needed for this index
                if tp_indx == 0:
                    model.c6.add( model.var_sats_estore[sat_indx,timepoint] ==  model.par_sats_estore_initial[sat_indx])
                    continue 

                #  minimum and maximum storage constraints
                model.c6.add( model.var_sats_estore[sat_indx,timepoint] >= model.par_sats_estore_min[sat_indx])
                model.c6.add( model.var_sats_estore[sat_indx,timepoint] <= model.par_sats_estore_max[sat_indx])

                # determine activity energy consumption
                activity_delta_e = 0 
                #  get the activities that were active during the time step immediately preceding time point
                activities = activity_dancecards[sat_indx].get_objects_pre_timepoint_indx(timepoint)
                for act in activities:
                    act_uindx = all_acts_by_obj[act]
                    activity_delta_e += (
                        model.par_sats_edot_by_act[sat_indx][self.act_type_map[type(act)]] 
                        * model.var_activity_utilization[act_uindx]
                        * model.par_resource_delta_t
                    )

                # determine whether or not the satellite is in sunlight and thus charging
                charging = True
                charging_delta_e = model.par_sats_edot_by_act[sat_indx]['charging']*model.par_resource_delta_t if charging else 0

                #  base-level satellite energy usage (not including additional activities)
                base_delta_e = model.par_sats_edot_by_act[sat_indx]['base']*model.par_resource_delta_t

                # maximum bound of energy at current time step based on last time step
                model.c6.add( model.var_sats_estore[sat_indx,timepoint] <= 
                    model.var_sats_estore[sat_indx,timepoint-1]
                    + activity_delta_e
                    + charging_delta_e
                    + base_delta_e
                )

                # maximum bound of energy at current time step based on last time step
                model.c6.add( model.var_sats_estore[sat_indx,timepoint] >= 
                    model.var_sats_estore[sat_indx,timepoint-1]
                    + activity_delta_e
                    + base_delta_e
                )





        ##############################
        #  Make objective
        ##############################

        def obj_rule(model):
            return (
                # model.par_obj_weight1 * 1/self.num_paths/model.par_obs_dv * sum(model.var_path_dv_dlnk[p,i,k] for i,k in model.dlnk_subscripts for p in model.paths )  +
                # model.par_obj_weight2 * 1/self.num_paths                  * sum(model.var_path_indic[p] for p in model.paths) +
                # model.par_obj_weight3 * 1/self.num_paths                  * sum(model.var_dlnk_path_occ[p,i,k]*model.par_dlnk_sf[i,k] for i,k in model.dlnk_subscripts for p in model.paths )
                sum(model.par_path_dv[p]*model.var_path_utilization[p] for p in model.paths) +
                sum(model.var_path_indic[p] for p in model.paths)
            )
        model.obj = pe.Objective( rule=obj_rule, sense=pe.maximize )

        self.model = model

    # taken in part from Jeff Menezes' code at https://github.mit.edu/jmenezes/Satellite-MILP/blob/master/sat_milp_pyomo.py
    def solve(self):

        solver = po.SolverFactory(self.solver_name)
        solver.options['timelimit'] = self.solver_max_runtime

        # if we're running things remotely, then we will use the NEOS server (https://neos-server.org/neos/)
        if self.solver_run_remotely:
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
        else:
            # something else is wrong
            print (results.solver)

    def print_sol_all(self):
        for v in self.model.component_objects(pe.Var, active=True):
            print ("Variable",v)
            varobject = getattr(self.model, str(v))
            for index in varobject:
                print (" ",index, varobject[index].value)

    def print_sol(self):
        for v in self.model.component_objects(pe.Var, active=True):
            if str (v) =='var_activity_utilization': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" ",index, val)

            elif str (v) =='var_path_utilization': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" ",index, val)
            
            elif str (v) =='var_path_indic': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" ",index, val)

    def extract_utilized_routes( self,verbose = False):
        #  note that routes are the same as paths

        # scheduled_dv
        scheduled_routes_flat = []

        if verbose:
            print ('utilized routes:')

        # figure out which paths were used, and add the scheduled data volume to each
        for p in self.model.paths:
            if pe.value(self.model.var_path_indic[p]) >= 1.0 - self.binary_epsilon:
                scheduled_route =  self.routes_flat[p]
                scheduled_route.scheduled_dv = pe.value(self.model.var_path_utilization[p])* self.model.par_path_dv[p]
                scheduled_routes_flat. append (scheduled_route)

                if verbose:
                    print(scheduled_route)

        return scheduled_routes_flat

    def extract_resource_usage( self, decimation_factor =1, verbose = False):
        #  note that routes are the same as paths

        energy_usage = {}

        t_vals = []
        e_vals = [[] for sat_indx in range ( self.num_sats)]

        # TODO: this code feels super inefficient somehow.  make it better?
        for sat_indx, sat in enumerate (self.model.sats):
            last_tp_indx = 0
            for tp_indx, tp in enumerate (self.model.timepoints):

                if (tp_indx - last_tp_indx) < decimation_factor:
                    continue
                else:
                    e_vals[sat_indx].append(pe.value(self.model.var_sats_estore[sat,tp]))

                    if sat_indx == 0:
                        t_vals.append(self.timepoints_m[tp])

        energy_usage['time_mins'] = t_vals
        energy_usage['e_sats'] = e_vals

        return  energy_usage


