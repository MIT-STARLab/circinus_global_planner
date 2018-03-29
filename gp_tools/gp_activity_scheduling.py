# contains model and Solver for global planner activity scheduling capability
# 
# @author Kit Kennedy
#
#  note that a path is the same as a route. 

from  datetime import timedelta
from copy import  deepcopy
from math import ceil

from pyomo import environ  as pe
from pyomo import opt  as po

from circinus_tools  import time_tools as tt
from .custom_activity_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow
from .routing_objects import DataRoute
from .schedule_objects import Dancecard

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

        scenario_params = gp_params['gp_orbit_prop_params']['scenario_params']
        sat_params = gp_params['gp_orbit_prop_params']['sat_params']
        as_params = gp_params['gp_general_params']['activity_scheduling_params']
        gp_inst_params = gp_params['gp_instance_params']['activity_scheduling_params']
        
        self.latency_params = gp_params['gp_general_params']['other_params']['latency_calculation']
        self.solver_max_runtime =as_params['solver_max_runtime_s']
        self.solver_name =as_params['solver_name']
        self.solver_optimality_gap =as_params['solver_optimality_gap']
        self.solver_run_remotely =as_params['solver_run_remotely']
        self.min_forked_route_dv =as_params['min_forked_route_dv_Mb']
        self.num_sats=sat_params['num_sats']
        self.transition_time_s=as_params['transition_time_s']
        self.sched_start_utc_dt  = gp_inst_params['start_utc_dt']
        self.sched_end_utc_dt  = gp_inst_params['end_utc_dt']
        self.resource_delta_t_s  =as_params['resource_delta_t_s']
        self.dv_epsilon = as_params['dv_epsilon_Mb']
        self.resource_margin_num_timepoints = as_params['resource_margin_num_timepoints']

        self.obj_weights =as_params['obj_weights']

        self.power_params = sat_params['power_params_sorted']
        self.initial_state = sat_params['initial_state_sorted']

        # these lists are in order of satellite index because we've sorted 
        self.sats_init_estate_Wh = [sat_state['batt_e_Wh'] for sat_state in self.initial_state]
        self.sats_emin_Wh = [p_params['battery_storage_Wh']['e_min'][p_params['battery_option']] for p_params in self.power_params]
        self.sats_emax_Wh = [p_params['battery_storage_Wh']['e_max'][p_params['battery_option']] for p_params in self.power_params]

        self.energy_unit = "Wh"

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
            XlnkWindow: 'xlnk',
            EclipseWindow: 'charging'
        }

        self.min_act_duration_s = {
            ObsWindow: as_params['min_duration_s']['obs'],
            DlnkWindow: as_params['min_duration_s']['dlnk'],
            XlnkWindow: as_params['min_duration_s']['xlnk']
        }

        self.standard_activities = [ObsWindow,DlnkWindow,XlnkWindow]

        if self.latency_params['obs'] not in ['start','end']:
            raise NotImplementedError
        if self.latency_params['dlnk'] not in ['start','end','center']:
            raise NotImplementedError

    def get_stats(self,verbose=True):
        stats = {}
        stats['num_routes'] = sum([len ( self.routes_flat)])
        stats['num_acts'] = sum([len ( self.all_acts_indcs)])
        stats['num_obs_acts'] = sum([len ( self.obs_act_indcs)])
        stats['num_link_acts'] = sum([len ( self.link_act_indcs)])
        stats['num_variables'] = self.model.nvariables ()
        stats['num_nobjectives'] = self.model.nobjectives ()
        stats['num_nconstraints'] = self.model.nconstraints ()

        if verbose:
            print ( "Considering %d routes" % (stats['num_routes']))
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
                    

    def filter_routes( self,routes_flat):

        new_routes = []
        for dr in routes_flat:
            dr_start = dr.get_obs().start
            dr_end = dr.get_dlnk().end

            if dr_start < self.sched_start_utc_dt or dr_end > self.sched_end_utc_dt:
                pass
            else:
                new_routes.append (dr)

        return new_routes

    @staticmethod
    def get_path_latency_score_factors(routes_flat,path_indcs_by_obs_act,latency_params):

        path_latency_sf_by_path_indx = {}

        # loop through all satellites and their downlinks
        # explicitly indexing by num_sats just to include a bit of error checking
        for obs,dr_indcs in path_indcs_by_obs_act.items():
            latencies = []
            for dr_indx in dr_indcs:
                dr = routes_flat[dr_indx]

                latencies.append(
                    dr.get_latency(
                        'minutes',
                        obs_option = latency_params['obs'], 
                        dlnk_option = latency_params['dlnk']
                    )
                )

             #  the shortest latency path for this observation has a score factor of 1.0, and the score factors for the other paths decrease as the inverse of increasing latency
            min_lat = min(latencies)
            for lat_indx, lat in enumerate(latencies):
                dr_indx = dr_indcs[lat_indx]
                path_latency_sf_by_path_indx[dr_indx] = min_lat/lat
 
        return path_latency_sf_by_path_indx

    def make_model ( self,routes_flat, ecl_winds, verbose = True):
        model = pe.ConcreteModel()

        # print(routes_flat)

        # filter the routes to make sure that  none of their activities fall outside the scheduling window
        routes_flat = self.filter_routes(routes_flat)

        self.routes_flat = routes_flat
        self.ecl_winds = ecl_winds

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

        path_latency_sf_by_path_indx =  self.get_path_latency_score_factors(
            routes_flat,
            path_indcs_by_obs_act,
            self.latency_params
        )

        # construct a set of dance cards for every satellite, 
        # each of which keeps track of all of the activities of satellite 
        # can possibly execute at any given time slice delta T. 
        act_dancecards = [Dancecard(self.sched_start_utc_dt,self.sched_end_utc_dt,self.resource_delta_t_s,mode='timestep') for sat_indx in range (self.num_sats)]
        for sat_indx in range (self.num_sats): 
            act_dancecards[sat_indx].add_winds_to_dancecard(sat_acts[sat_indx])
            act_dancecards[sat_indx].add_winds_to_dancecard(ecl_winds[sat_indx])

        self.all_acts_indcs = all_acts_indcs
        self.obs_act_indcs = obs_act_indcs
        self.link_act_indcs = link_act_indcs
        self.all_acts_by_obj = all_acts_by_obj

        #  subscript for each path p
        model.paths = pe.Set(initialize= range(len(routes_flat)))
        #  subscript for each activity a
        model.acts = pe.Set(initialize= all_acts_indcs)
        #  subscript for each satellite
        model.sats = pe.Set(initialize=  range ( self.num_sats))

        # timepoints is the indices, which starts at 0 
        #  NOTE: we assume the same time system for every satellite
        self.time_getter_dc = act_dancecards[0]
        num_timepoints = act_dancecards[0].num_timepoints
        model.timepoint_indcs = pe.Set(initialize=  self.time_getter_dc.get_tp_indcs())

        #  unique indices for observation and link acts
        model.obs_acts = pe.Set(initialize= obs_act_indcs)
        model.link_acts = pe.Set(initialize= link_act_indcs)

        ##############################
        #  Make parameters
        ##############################

        model.par_min_forked_route_dv = pe.Param (initialize=self.min_forked_route_dv)
        model.par_obs_dv = pe.Param(model.obs_acts,initialize =dv_by_obs_act)
        model.par_total_obs_dv = sum(dv_by_obs_act.values())
        model.par_link_dv = pe.Param(model.link_acts,initialize =dv_by_link_act)
        model.par_act_dv = pe.Param(model.acts,initialize =dv_by_act)
        model.par_path_dv = pe.Param(model.paths,initialize ={ i: dr.data_vol for i,dr in enumerate (routes_flat)})
        # each of these is essentially a dictionary indexed by link or obs act indx, with  the value being a list of path indices that are included within that act
        # these are valid indices into model.paths
        model.par_path_subscrs_by_link_act = pe.Param(model.link_acts,initialize =path_indcs_by_link_act)
        model.par_path_subscrs_by_obs_act = pe.Param(model.obs_acts,initialize =path_indcs_by_obs_act)
        model.par_path_subscrs_by_act = pe.Param(model.acts,initialize =path_indcs_by_act)


        if self.energy_unit == "Wh":
            model.par_resource_delta_t = pe.Param (initialize= self.resource_delta_t_s/3600)
        else: 
            raise NotImplementedError
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
        model.var_sats_estore  = pe.Var (model.sats,  model.timepoint_indcs,  within = pe.NonNegativeReals)
        

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
            return model.par_path_dv[p]*model.var_path_utilization[p] >= model.par_min_forked_route_dv*model.var_path_indic[p]
        model.c2 =pe.Constraint ( model.paths,  rule=c2_rule)

        def c3_rule( model,a):
            return model.var_act_indic[a] >=  model.var_activity_utilization[a]
        model.c3 =pe.Constraint ( model.acts,  rule=c3_rule)  

        def c3c_rule( model,p):
            return model.var_path_indic[p] >=  model.var_path_utilization[p]
        model.c3c =pe.Constraint ( model.paths,  rule=c3c_rule)

        #  activity overlap constraints [4],[5],[5b]
        model.c4  = pe.ConstraintList()
        model.c5  = pe.ConstraintList()
        model.c5b  = pe.ConstraintList()
        for sat_indx in range (self.num_sats):
            num_sat_acts = len(sat_acts[sat_indx])
            for  first_act_indx in  range (num_sat_acts):

                act1 = sat_acts[sat_indx][first_act_indx]
                # get the unique index into model.acts
                act1_uindx = all_acts_by_obj[act1]
                length_1 = model.par_act_dv[act1_uindx]*model.var_activity_utilization[act1_uindx]/act1.ave_data_rate
                model.c5b.add( length_1 >= model.var_act_indic[act1_uindx] * self.min_act_duration_s[type(act1)])
                
                for  second_act_indx in  range (first_act_indx+1,num_sat_acts):
                    act2 = sat_acts[sat_indx][ second_act_indx]
                    # get the unique index into model.acts
                    act2_uindx = all_acts_by_obj[act2]

                    # if there is enough transition time between the two activities, no constraint needs to be added
                    if (act2.start - act1.end).total_seconds() >= self.transition_time_s['simple']:
                        #  don't need to do anything,  continue on to next activity pair
                        continue

                    #  if the activities overlap in center time, then it's not possible to have sufficient transition time between them
                    #  add constraint to rule out the possibility of scheduling both of them
                    elif (act2.center - act1.center).total_seconds() <= self.transition_time_s['simple']:
                        model.c4.add( model.var_act_indic[act1_uindx]+ model.var_act_indic[act2_uindx] <= 1)

                    # If they don't overlap in center time, but they do overlap to some amount, then we need to constrain their end and start times to be consistent with one another
                    else:
                        center_time_diff = (act2.center - act1.center).total_seconds()
                        # this is the adjustment added to the center time to get to the start or end of the activity
                        time_adjust_1 = model.par_act_dv[act1_uindx]*model.var_activity_utilization[act1_uindx]/2/act1.ave_data_rate
                        time_adjust_2 = model.par_act_dv[act2_uindx]*model.var_activity_utilization[act2_uindx]/2/act2.ave_data_rate
                        model.c5.add( center_time_diff - time_adjust_1 - time_adjust_2 >= self.transition_time_s['simple'])

        #  energy constraints [6]
        model.c6  = pe.ConstraintList()
        for sat_indx in range (self.num_sats): 

            # tp_indx serves as an index into the satellite activity dance cards
            for tp_indx in model.timepoint_indcs:
                #  constraining first time step to initial energy storage
                #  continue for loop afterwards because no geq/leq constraints needed for this index
                if tp_indx == 0:
                    model.c6.add( model.var_sats_estore[sat_indx,tp_indx] ==  model.par_sats_estore_initial[sat_indx])
                    continue 

                #  minimum and maximum storage constraints
                model.c6.add( model.var_sats_estore[sat_indx,tp_indx] >= model.par_sats_estore_min[sat_indx])
                model.c6.add( model.var_sats_estore[sat_indx,tp_indx] <= model.par_sats_estore_max[sat_indx])

                # determine activity energy consumption
                charging = True
                activity_delta_e = 0 
                #  get the activities that were active during the time step immediately preceding time point
                activities = act_dancecards[sat_indx].get_objects_at_ts_pre_tp_indx(tp_indx)
                for act in activities:
                    #  if this is a "standard activity" that we can choose to perform or not
                    if type(act) in self.standard_activities:
                        act_uindx = all_acts_by_obj[act]
                        activity_delta_e += (
                            model.par_sats_edot_by_act[sat_indx][self.act_type_map[type(act)]] 
                            * model.var_activity_utilization[act_uindx]
                            * model.par_resource_delta_t
                        )

                    #  if the satellite is not in sunlight then we can't charge
                    elif type(act) == EclipseWindow:
                        charging = False

                # add in charging energy contribution ( if possible)
                charging_delta_e = model.par_sats_edot_by_act[sat_indx]['charging']*model.par_resource_delta_t if charging else 0

                #  base-level satellite energy usage (not including additional activities)
                base_delta_e = model.par_sats_edot_by_act[sat_indx]['base']*model.par_resource_delta_t

                # maximum bound of energy at current time step based on last time step
                model.c6.add( model.var_sats_estore[sat_indx,tp_indx] <= 
                    model.var_sats_estore[sat_indx,tp_indx-1]
                    + activity_delta_e
                    + charging_delta_e
                    + base_delta_e
                )

                # maximum bound of energy at current time step based on last time step
                model.c6.add( model.var_sats_estore[sat_indx,tp_indx] >= 
                    model.var_sats_estore[sat_indx,tp_indx-1]
                    + activity_delta_e
                    + base_delta_e
                )

                #  TODO: reassess if adding this version of the code back in is faster
                # if charging:
                #     # maximum bound of energy at current time step based on last time step
                #     model.c6.add( model.var_sats_estore[sat_indx,tp_indx] <= 
                #         model.var_sats_estore[sat_indx,tp_indx-1]
                #         + activity_delta_e
                #         + charging_delta_e
                #         + base_delta_e
                #     )

                #     # maximum bound of energy at current time step based on last time step
                #     model.c6.add( model.var_sats_estore[sat_indx,tp_indx] >= 
                #         model.var_sats_estore[sat_indx,tp_indx-1]
                #         + activity_delta_e
                #         + base_delta_e
                #     )
                # else:
                #     # energy at time step is deterministic function of last time step
                #     model.c6.add( model.var_sats_estore[sat_indx,tp_indx] == 
                #         model.var_sats_estore[sat_indx,tp_indx-1]
                #         + activity_delta_e
                #         + base_delta_e
                #     )






        ##############################
        #  Make objective
        ##############################

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        #  determine which time points to use as "spot checks" on resource margin. These are the points that will be used in the objective function for maximizing resource margin
        timepoint_spacing = ceil(num_timepoints/self.resource_margin_num_timepoints)
        # need to turn the generator into a list for slicing
        #  note: have to get the generator again
        decimated_tp_indcs = list(self.time_getter_dc.get_tp_indcs())[::timepoint_spacing]
        rsrc_norm_f = len(decimated_tp_indcs) * len(model.sats)

        def obj_rule(model):
            return (
                self.obj_weights['obs_dv'] * 1/model.par_total_obs_dv * sum(model.par_path_dv[p]*model.var_path_utilization[p] for p in model.paths) 
                + self.obj_weights['route_latency'] * 1/len(model.obs_acts) * sum(path_latency_sf_by_path_indx[p]*model.var_path_utilization[p] for p in model.paths)
                + self.obj_weights['energy_storage'] * 1/rsrc_norm_f * sum(model.var_sats_estore[sat_indx,tp_indx]/model.par_sats_estore_max[sat_indx] for tp_indx in decimated_tp_indcs for sat_indx in model.sats)
            )
        model.obj = pe.Objective( rule=obj_rule, sense=pe.maximize )

        self.model = model

    # taken in part from Jeff Menezes' code at https://github.mit.edu/jmenezes/Satellite-MILP/blob/master/sat_milp_pyomo.py
    def solve(self):

        solver = po.SolverFactory(self.solver_name)
        solver.options['timelimit'] = self.solver_max_runtime
        if self.solver_name == 'gurobi':
            # note default for this is 1e-4, or 0.01%
            solver.options['MIPGap'] = self.solver_optimality_gap

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

    def extract_utilized_routes( self, copy_windows = True, verbose = False):
        #  note that routes are the same as paths

        # scheduled_dv
        scheduled_routes_flat = []

        if verbose:
            print ('utilized routes:')

        def copy_choice(wind):
            if copy_windows:
                return deepcopy(wind)
            else:
                return wind

        # figure out which paths were used, and add the scheduled data volume to each
        for p in self.model.paths:
            if pe.value(self.model.var_path_indic[p]) >= 1.0 - self.binary_epsilon:
                scheduled_route =  copy_choice (self.routes_flat[p]) 
                scheduled_route.scheduled_dv = pe.value(self.model.var_path_utilization[p])* self.model.par_path_dv[p]
                scheduled_routes_flat. append (scheduled_route)

                if verbose:
                    print(scheduled_route)

        # examine the schedulable data volume for every activity window, checking as we go that the data volume is sufficient for at least the route in which the window is found
        #  note that this code is slightly inefficient because it might duplicate windows across routes. that's fine though, because we're thorough in checking across all routes
        for dr in scheduled_routes_flat:
            for wind in dr.route:
                act_indx = self.all_acts_by_obj[wind]
                w_s_dv = wind.data_vol * pe.value(self.model.var_act_indic[act_indx])
                #  initialize this while we're here
                wind.scheduled_data_vol = 0
                if w_s_dv < dr.scheduled_dv - self.dv_epsilon: 
                    raise RuntimeWarning('inconsistent activity scheduling results: activity data volume (%f) smaller than route data volume (%f) [%s,%s]'%(w_s_dv,dr.scheduled_dv,dr,wind))


        #  now we want to mark the real scheduled data volume for every window. We need to do this separately because the model.var_act_indic continuous variables only give an upper bound on the data volume for an activity. we only actually need to use as much data volume as the data routes want to push through the window
        #  add data volume for every route passing through every window
        for dr in scheduled_routes_flat:
            for wind in dr.route:
                wind.scheduled_data_vol += dr.scheduled_dv

        # update the window beginning and end times based upon their amount of scheduled data volume
        # keep track of which ones we've updated, because we should only update once
        updated_winds = set()
        for dr in scheduled_routes_flat:
            for wind in dr.route:
                if not wind in updated_winds:
                    # note that the line below seems like it may break the scheduled times for activities by specifying a minimum activity duration. however, this minimum activity duration is already accounted for in scheduling constraints
                    wind.update_duration_from_scheduled_dv (min_duration_s=self.min_act_duration_s[type(wind)])
                    updated_winds.add(wind)

        return scheduled_routes_flat

    def extract_resource_usage( self, decimation_factor =1, verbose = False):
        #  note that routes are the same as paths

        energy_usage = {}

        t_vals = []
        e_vals = [[] for sat_indx in range ( self.num_sats)]

        # TODO: this code feels super inefficient somehow.  make it better?
        for sat_indx, sat in enumerate (self.model.sats):
            last_tp_indx = 0
            for tp_indx in self.model.timepoint_indcs:

                if (tp_indx - last_tp_indx) < decimation_factor:
                    continue
                else:
                    e_vals[sat_indx].append(pe.value(self.model.var_sats_estore[sat,tp_indx]))

                    if sat_indx == 0:
                        t_vals.append(self.time_getter_dc.get_tp_from_tp_indx(tp_indx,out_units='minutes'))

        energy_usage['time_mins'] = t_vals
        energy_usage['e_sats'] = e_vals

        return  energy_usage


