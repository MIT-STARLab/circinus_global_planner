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
from .custom_activity_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow
from .routing_objects import DataMultiRoute
from .schedule_objects import Dancecard
from circinus_tools  import  constants as const

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
        self.solver_name =as_params['solver_name']
        self.solver_params =as_params['solver_params']
        self.min_as_route_dv =as_params['min_as_route_dv_Mb']
        self.num_sats=sat_params['num_sats']
        self.sched_start_utc_dt  = gp_inst_params['start_utc_dt']
        self.sched_end_utc_dt  = gp_inst_params['end_utc_dt']
        self.resource_delta_t_s  =as_params['resource_delta_t_s']
        self.enforce_energy_storage_constr  =as_params['enforce_energy_storage_constr']
        self.enforce_data_storage_constr  =as_params['enforce_data_storage_constr']

        #  the "effectively zero" number.
        self.dv_epsilon = as_params['dv_epsilon_Mb']
        self.resource_margin_obj_num_timepoints = as_params['resource_margin_obj_num_timepoints']

        self.obj_weights =as_params['obj_weights']

        self.use_symmetric_xlnk_windows = gp_params['gp_general_params']['other_params']['use_symmetric_xlnk_windows']

        self.power_params = sat_params['power_params_sorted']
        self.data_storage_params = sat_params['data_storage_params_sorted']
        self.initial_state = sat_params['initial_state_sorted']
        sat_activity_params = sat_params['activity_params']

        # these lists are in order of satellite index because we've sorted 
        self.sats_init_estate_Wh = [sat_state['batt_e_Wh'] for sat_state in self.initial_state]
        self.sats_emin_Wh = [p_params['battery_storage_Wh']['e_min'][p_params['battery_option']] for p_params in self.power_params]
        self.sats_emax_Wh = [p_params['battery_storage_Wh']['e_max'][p_params['battery_option']] for p_params in self.power_params]

        self.sats_dmin_Mb = [1000*ds_params['data_storage_Gbit']['d_min'][ds_params['storage_option']] for ds_params in self.data_storage_params]
        self.sats_dmax_Mb = [1000*ds_params['data_storage_Gbit']['d_max'][ds_params['storage_option']] for ds_params in self.data_storage_params]

        self.energy_unit = "Wh"

        self.sats_edot_by_act_W = []
        for p_params in self.power_params:
            sat_edot_by_act = {}
            sat_edot_by_act['base'] = p_params['power_consumption_W']['base'][p_params['base_option']]
            sat_edot_by_act['obs'] = p_params['power_consumption_W']['obs'][p_params['obs_option']]
            sat_edot_by_act['dlnk'] = p_params['power_consumption_W']['dlnk'][p_params['dlnk_option']]
            sat_edot_by_act['xlnk-tx'] = p_params['power_consumption_W']['xlnk-tx'][p_params['xlnk_tx_option']]
            sat_edot_by_act['xlnk-rx'] = p_params['power_consumption_W']['xlnk-rx'][p_params['xlnk_rx_option']]
            sat_edot_by_act['charging'] = p_params['power_consumption_W']['orbit_insunlight_average_charging'][p_params['charging_option']]
            self.sats_edot_by_act_W.append (sat_edot_by_act)

        self.act_type_map = {
            ObsWindow: 'obs',
            DlnkWindow: 'dlnk',
            'xlnk-tx': 'xlnk-tx',
            'xlnk-rx': 'xlnk-rx',
            EclipseWindow: 'charging'
        }

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
        self.standard_activities = [ObsWindow,DlnkWindow]

        #  this should be as small as possible to prevent ill conditioning, but big enough that score factor constraints are still valid. 
        #  note: the size of this value is checked in make_model() below
        self.big_M_lat = 1e6


        # this big M should be conformatably bigger than the duration of any activity. It's used to switch on/off constraints that are related to differences in start/end times of activities, and in the worst case we expect such differences to max out at the legnth of the longest-duration activity. As long as we provide some comfortable margin past that duration, should be fine.
        # setting to 20 minutes
        self.big_M_act_t_dur_s = 30*60

        # this big M should be twice as large as any activity data volume 
        self.big_M_act_dv = 200000 # in Mb

    def get_stats(self,verbose=True):

        num_winds_per_route = [len(dmr.get_winds()) for dmr in self.routes_flat]
        num_routes_by_act = {act:len(self.dmr_indcs_by_act[act_indx]) for act,act_indx in self.all_acts_by_obj.items()}

        stats = {}
        stats['num_routes'] = sum([len ( self.routes_flat)])
        stats['num_acts'] = sum([len ( self.all_acts_indcs)])
        stats['num_obs_acts'] = sum([len ( self.obs_act_indcs)])
        stats['num_link_acts'] = sum([len ( self.link_act_indcs)])
        stats['ave_num_winds_per_route'] = np.mean(num_winds_per_route)
        stats['max_num_winds_per_route'] = np.max(num_winds_per_route)
        stats['med_num_winds_per_route'] = np.median(num_winds_per_route)
        stats['num_variables'] = self.model.nvariables ()
        stats['num_nobjectives'] = self.model.nobjectives ()
        stats['num_nconstraints'] = self.model.nconstraints ()
        stats['num_dlnks'] = sum([1 for act in num_routes_by_act.keys() if type(act) == DlnkWindow])
        stats['num_xlnks'] = sum([1 for act in num_routes_by_act.keys() if type(act) == XlnkWindow])

        stats['ave_num_routes_per_act'] = np.mean(list(num_routes_by_act.values()))
        stats['ave_num_routes_per_obs'] = np.mean([num for act,num in num_routes_by_act.items() if type(act) == ObsWindow])
        stats['ave_num_routes_per_dlnk'] = np.mean([num for act,num in num_routes_by_act.items() if type(act) == DlnkWindow])
        stats['max_num_routes_per_dlnk'] = np.max([num for act,num in num_routes_by_act.items() if type(act) == DlnkWindow])
        stats['min_num_routes_per_dlnk'] = np.min([num for act,num in num_routes_by_act.items() if type(act) == DlnkWindow])
        stats['med_num_routes_per_dlnk'] = np.median([num for act,num in num_routes_by_act.items() if type(act) == DlnkWindow])
        stats['ave_num_routes_per_xlnk'] = np.mean([num for act,num in num_routes_by_act.items() if type(act) == XlnkWindow])
        stats['max_num_routes_per_xlnk'] = np.max([num for act,num in num_routes_by_act.items() if type(act) == XlnkWindow])


        if verbose:
            print ( "Considering %d routes" % (stats['num_routes']))
            print ( "Considering %d activity windows" % (stats['num_acts']))
            print ( "Considering %d observation windows" % (stats['num_obs_acts']))
            print ( "Considering %d link windows" % (stats['num_link_acts']))
            print ( "Considering %d xlink windows" % (stats['num_xlnks']))
            print ( "Considering %d dlink windows" % (stats['num_dlnks']))
            print ( "ave_num_winds_per_route: %0.2f"%(stats['ave_num_winds_per_route']))
            print ( "max_num_winds_per_route: %0.2f"%(stats['max_num_winds_per_route']))
            print ( "med_num_winds_per_route: %0.2f"%(stats['med_num_winds_per_route']))
            print ( "ave_num_routes_per_act: %0.2f"%(stats['ave_num_routes_per_act']))
            print ( "ave_num_routes_per_obs: %0.2f"%(stats['ave_num_routes_per_obs']))
            print ( "ave_num_routes_per_dlnk: %0.2f"%(stats['ave_num_routes_per_dlnk']))
            print ( "max_num_routes_per_dlnk: %0.2f"%(stats['max_num_routes_per_dlnk']))
            print ( "min_num_routes_per_dlnk: %0.2f"%(stats['min_num_routes_per_dlnk']))
            print ( "med_num_routes_per_dlnk: %0.2f"%(stats['med_num_routes_per_dlnk']))
            print ( "ave_num_routes_per_xlnk: %0.2f"%(stats['ave_num_routes_per_xlnk']))
            print ( "max_num_routes_per_xlnk: %0.2f"%(stats['max_num_routes_per_xlnk']))
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
        # dmr is data multi-route
        dmr_indcs_by_link_act = {}
        dmr_indcs_by_obs_act = {}
        dmr_indcs_by_act = {}
        dv_by_link_act = {}
        dv_by_obs_act = {}
        dv_by_act = {}

        sat_acts = [[] for sat_indx in range (self.num_sats)]
        sat_dlnks = [[] for sat_indx in range (self.num_sats)]

        new_act_indx = 0
        for dmr_indx, dmr in enumerate (routes_flat):
            for act in dmr.get_winds():

                # if we haven't yet seen this activity, then add it to bookkeeping
                if not act in all_acts_by_obj.keys():
                    act_indx = new_act_indx
                    all_acts_indcs.append(act_indx)
                    all_acts_by_indx[act_indx] = act
                    all_acts_by_obj[act] = act_indx
                    dmr_indcs_by_act[act_indx] = []
                    dmr_indcs_by_act[act_indx].append (dmr_indx)
                    dv_by_act[act_indx] = act.data_vol
                    new_act_indx += 1

                    # also need to add it to the list and dictionary for observations
                    if type(act) == ObsWindow:
                        sat_acts[act.sat_indx].append(act)
                        obs_act_indcs.append(act_indx)
                        dmr_indcs_by_obs_act[act_indx] = []
                        dmr_indcs_by_obs_act[act_indx].append (dmr_indx)
                        dv_by_obs_act[act_indx] = act.data_vol

                    # also need to add it to the list and dictionary for links
                    if type(act) == DlnkWindow:
                        link_act_indcs.append(act_indx)
                        dmr_indcs_by_link_act[act_indx] = []
                        dmr_indcs_by_link_act[act_indx].append (dmr_indx)
                        dv_by_link_act[act_indx] = act.data_vol
                        sat_acts[act.sat_indx].append(act)
                        # grab the dlnks for each sat too, while we're looping through
                        sat_dlnks[act.sat_indx].append(act)

                    if type(act) == XlnkWindow:
                        link_act_indcs.append(act_indx)
                        dmr_indcs_by_link_act[act_indx] = []
                        dmr_indcs_by_link_act[act_indx].append (dmr_indx)
                        dv_by_link_act[act_indx] = act.data_vol
                        sat_acts[act.sat_indx].append(act)
                        sat_acts[act.xsat_indx].append(act)

                #  if we have already seen the activity,  then just need to update the appropriate structures
                else:
                    act_indx = all_acts_by_obj[act]
                    dmr_indcs_by_act[act_indx].append (dmr_indx)

                    # add the data route index
                    if type(act) == ObsWindow:
                        dmr_indcs_by_obs_act[act_indx].append (dmr_indx)
                    if type(act) == DlnkWindow or type(act) == XlnkWindow:
                        dmr_indcs_by_link_act[act_indx].append (dmr_indx)

        # print (dmr_indcs_by_link_act)
        # print (dmr_indcs_by_obs_act)

        #  sort the activities, because we'll need that for constructing constraints
        for sat_indx in range (self.num_sats):
            sat_acts[sat_indx].sort(key=lambda x: x.center)
            sat_dlnks[sat_indx].sort(key=lambda x: x.center)

        return sat_acts,sat_dlnks,all_acts_indcs,dmr_indcs_by_act,dv_by_act,all_acts_by_obj,all_acts_by_indx,obs_act_indcs,dmr_indcs_by_obs_act,dv_by_obs_act,link_act_indcs,dmr_indcs_by_link_act,dv_by_link_act
                    

    def filter_routes( self,routes_flat):
        """ todo: should filter routes based on having the observation be within the window of interest for scheduling"""

        new_routes = []
        for dmr in routes_flat:
            dmr_start = dmr.get_obs().start
            dmr_end = dmr.get_dlnk().end

            if dmr_start < self.sched_start_utc_dt or dmr_end > self.sched_end_utc_dt:
                pass
            if dmr.get_dlnk().duration.total_seconds() < self.min_act_duration_s[DlnkWindow]:
                print('discarding too short dlnk window')
                pass
            else:
                new_routes.append (dmr)

        return new_routes

    @staticmethod
    def get_dmr_latency_score_factors(routes_flat,dmr_indcs_by_obs_act,latency_params):

        dmr_latency_sf_by_dmr_indx = {}

        # loop through all satellites and their downlinks
        # explicitly indexing by num_sats just to include a bit of error checking
        for obs,dmr_indcs in dmr_indcs_by_obs_act.items():
            latencies = []
            for dmr_indx in dmr_indcs:
                dmr = routes_flat[dmr_indx]

                latencies.append(
                    dmr.get_latency(
                        'minutes',
                        obs_option = latency_params['obs'], 
                        dlnk_option = latency_params['dlnk']
                    )
                )

             #  the shortest latency dmr (DataMultiRoute) for this observation has a score factor of 1.0, and the score factors for the other dmrs decrease as the inverse of increasing latency
            min_lat = min(latencies)
            for lat_indx, lat in enumerate(latencies):
                dmr_indx = dmr_indcs[lat_indx]
                dmr_latency_sf_by_dmr_indx[dmr_indx] = min_lat/lat
 
        return dmr_latency_sf_by_dmr_indx

    def make_model ( self,routes_flat, ecl_winds, verbose = True):
        model = pe.ConcreteModel()

        # filter the routes to make sure that  none of their activities fall outside the scheduling window
        routes_flat = self.filter_routes(routes_flat)

        self.routes_flat = routes_flat
        self.ecl_winds = ecl_winds

        #  should only be using data multi-route objects for activity scheduling, even if they're just a shallow wrapper around a DataRoute
        for dmr in routes_flat:
            assert(type(dmr) == DataMultiRoute)

        #  really useful code below!!!
        if verbose:
            pass

        ##############################
        #  Make indices/ subscripts
        ##############################

        try:
            (sat_acts,
                sat_dlnks,
                all_acts_indcs,
                dmr_indcs_by_act,
                dv_by_act,
                all_acts_by_obj,
                all_acts_by_indx,
                obs_act_indcs,
                dmr_indcs_by_obs_act,
                dv_by_obs_act,
                link_act_indcs,
                dmr_indcs_by_link_act,
                dv_by_link_act) =  self.get_activity_structs(routes_flat)

            dmr_latency_sf_by_dmr_indx =  self.get_dmr_latency_score_factors(
                routes_flat,
                dmr_indcs_by_obs_act,
                self.latency_params
            )

            # construct a set of dance cards for every satellite, 
            # each of which keeps track of all of the activities of satellite 
            # can possibly execute at any given time slice delta T. 
            # this is for constructing energy storage constraints
            # using resource_delta_t_s because this dancecard is solely for use in constructing resource constraints
            es_act_dancecards = [Dancecard(self.sched_start_utc_dt,self.sched_end_utc_dt,self.resource_delta_t_s,item_init=None,mode='timestep') for sat_indx in range (self.num_sats)]
            
            for sat_indx in range (self.num_sats): 
                es_act_dancecards[sat_indx].add_winds_to_dancecard(sat_acts[sat_indx])
                es_act_dancecards[sat_indx].add_winds_to_dancecard(ecl_winds[sat_indx])

            # this is for data storage
            # for each sat/timepoint, we store a list of all those data multi routes that are storing data on the sat at that timepoint
            ds_route_dancecards = [Dancecard(self.sched_start_utc_dt,self.sched_end_utc_dt,self.resource_delta_t_s,item_init=None,mode='timepoint') for sat_indx in range (self.num_sats)]
            
            # add data routes to the dancecard
            for dmr_indx,dmr in enumerate(routes_flat):
                # list of type routing_objects.SatStorageInterval
                dmr_ds_intervs = dmr.get_data_storage_intervals()

                for interv in dmr_ds_intervs:
                    # store the dmr object at this timepoint
                    ds_route_dancecards[interv.sat_indx].add_item_in_interval(dmr_indx,interv.start,interv.end)

        except IndexError:
            raise RuntimeWarning('sat_indx out of range. Are you sure all of your input files are consistent? (including pickles)')        
        self.dmr_indcs_by_act = dmr_indcs_by_act
        self.all_acts_indcs = all_acts_indcs
        self.obs_act_indcs = obs_act_indcs
        self.link_act_indcs = link_act_indcs
        self.all_acts_by_obj = all_acts_by_obj
        self.all_acts_by_indx = all_acts_by_indx

        # these subscripts probably should've been done using the unique IDs for the objects, rather than their arbitrary locations within a list. Live and learn, hélas...

        #  subscript for each dmr (data multi route) p  (p index is a hold-over from when I called them paths)
        model.dmrs = pe.Set(initialize= range(len(routes_flat)))
        #  subscript for each activity a
        model.acts = pe.Set(initialize= all_acts_indcs)
        #  subscript for each satellite
        model.sats = pe.Set(initialize=  range ( self.num_sats))

        # timepoints is the indices, which starts at 0 
        #  NOTE: we assume the same time system for every satellite
        self.es_time_getter_dc = es_act_dancecards[0]
        es_num_timepoints = es_act_dancecards[0].num_timepoints
        model.es_timepoint_indcs = pe.Set(initialize=  self.es_time_getter_dc.get_tp_indcs())

        self.ds_time_getter_dc = ds_route_dancecards[0]
        model.ds_timepoint_indcs = pe.Set(initialize=  self.ds_time_getter_dc.get_tp_indcs())

        #  unique indices for observation and link acts
        model.obs_acts = pe.Set(initialize= obs_act_indcs)
        model.link_acts = pe.Set(initialize= link_act_indcs)

        if self.solver_name == 'gurobi' or self.solver_name == 'cplex':
            int_feas_tol = self.solver_params['integer_feasibility_tolerance']
        else:
            raise NotImplementedError

        for p,lat_sf in dmr_latency_sf_by_dmr_indx.items():        
            if lat_sf > int_feas_tol*self.big_M_lat:
                raise RuntimeWarning('big_M_lat (%f) is not large enough for latency score factor %f and integer feasibility tolerance %f (dmr index %d)'%(self.big_M_lat,lat_sf,int_feas_tol,p))

        for act_obj in all_acts_by_obj.keys():
            if 2*(act_obj.end-act_obj.start).total_seconds() > self.big_M_act_t_dur_s:
                raise RuntimeWarning('big_M_act_t_dur_s (%f) is not large enough for act of duration %s and integer feasibility tolerance %f (act string %s)'%(self.big_M_act_t_dur_s,act_obj.end-act_obj.start,int_feas_tol,act_obj))
            # if 2*(act_obj.end-act_obj.start).total_seconds() > (1-int_feas_tol) *  self.big_M_act_t_dur_s:
                # raise RuntimeWarning('big_M_act_t_dur_s (%f) is not large enough for act of duration %s and integer feasibility tolerance %f (act string %s)'%(self.big_M_act_t_dur_s,act_obj.end-act_obj.start,int_feas_tol,act_obj))

            if 2*act_obj.data_vol > (1-int_feas_tol) * self.big_M_act_dv:
                raise RuntimeWarning('big_M_act_dv (%f) is not large enough for act of dv %f and integer feasibility tolerance %f (act string %s)'%(self.big_M_act_dv,act_obj.data_vol,int_feas_tol,act_obj))
                


        ##############################
        #  Make parameters
        ##############################

        model.par_min_as_route_dv = pe.Param (initialize=self.min_as_route_dv)
        model.par_obs_dv = pe.Param(model.obs_acts,initialize =dv_by_obs_act)
        model.par_total_obs_dv = sum(dv_by_obs_act.values())
        model.par_link_dv = pe.Param(model.link_acts,initialize =dv_by_link_act)
        model.par_act_dv = pe.Param(model.acts,initialize =dv_by_act)
        #  data volume for each data multi-route
        model.par_dmr_dv = pe.Param(model.dmrs,initialize ={ dmr_indx: dmr.data_vol for dmr_indx,dmr in enumerate (routes_flat)})
        #  data volume for each activity in each data multi-route

        # todo: will probably have to fix this to not be the full set multiplication of model.dmrs, model.acts - can use dmr_indcs_by_act
        model.par_dmr_act_dv = pe.Param(
            model.dmrs,
            model.acts,
            initialize = { (dmr_indx,self.all_acts_by_obj[act]): 
                dmr.data_vol_for_wind(act) for dmr_indx,dmr in enumerate (routes_flat) for act in dmr.get_winds()
            }
        )

        # each of these is essentially a dictionary indexed by link or obs act indx, with  the value being a list of dmr indices that are included within that act
        # these are valid indices into model.dmrs
        model.par_dmr_subscrs_by_link_act = pe.Param(model.link_acts,initialize =dmr_indcs_by_link_act)
        model.par_dmr_subscrs_by_obs_act = pe.Param(model.obs_acts,initialize =dmr_indcs_by_obs_act)
        model.par_dmr_subscrs_by_act = pe.Param(model.acts,initialize =dmr_indcs_by_act)


        if self.energy_unit == "Wh":
            model.par_resource_delta_t = pe.Param (initialize= self.resource_delta_t_s/3600)
        else: 
            raise NotImplementedError
        model.par_sats_estore_initial = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_init_estate_Wh)})
        model.par_sats_estore_min = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_emin_Wh)})
        model.par_sats_estore_max = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_emax_Wh)})
        model.par_sats_edot_by_act = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_edot_by_act_W)})

        model.par_sats_dstore_min = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_dmin_Mb)})
        model.par_sats_dstore_max = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_dmax_Mb)})

        ##############################
        #  Make variables
        ##############################


        # activity utilization variable indicating how much of an activity's capacity is used [1]
        model.var_activity_utilization  = pe.Var (model.acts, bounds =(0,1))
        # dmr utilization variable indicating how much of a dmr's capacity is used [2]
        model.var_dmr_utilization  = pe.Var (model.dmrs, bounds =(0,1))
        #  indicator variables for whether or not dmrs [3] and activities [4] have been chosen
        model.var_dmr_indic  = pe.Var (model.dmrs, within = pe.Binary)
        # model.var_dmr_indic  = pe.Var (model.dmrs, bounds =(0,1))
        model.var_act_indic  = pe.Var (model.acts, within = pe.Binary)
        # model.var_act_indic  = pe.Var (model.acts, bounds =(0,1))
        
        # satellite energy storage
        model.var_sats_estore  = pe.Var (model.sats,  model.es_timepoint_indcs,  within = pe.NonNegativeReals)

        # satellite data storage (data buffers)
        model.var_sats_dstore  = pe.Var (model.sats,  model.ds_timepoint_indcs,  within = pe.NonNegativeReals)

        model.var_dmr_latency_sf_by_obs_indx = pe.Var (model.obs_acts,  within = pe.NonNegativeReals)
        

        ##############################
        #  Make constraints
        ##############################

        # TODO: renumber  these with the final numbering

        # note that the observations show up within model.acts as well, so we also constraint route scheduled DV by the real available DV from each observation
        def c1_rule( model,a):
            return (model.par_act_dv[a]*model.var_activity_utilization[a] -
                        sum(model.par_dmr_act_dv[p,a]*model.var_dmr_utilization[p] 
                            for p in model.par_dmr_subscrs_by_act[a]) 
                    >= 0)
        model.c1 =pe.Constraint ( model.acts,  rule=c1_rule)

        model.c1b  = pe.ConstraintList()
        for a in model.acts:
            for p in model.par_dmr_subscrs_by_act[a]:
                model.c1b.add(model.var_act_indic[a] >= model.var_dmr_indic[p]) 

        model.c1c  = pe.ConstraintList()
        for act in sat_acts:
            if type(act) == DlnkWindow:
                act_uindx = all_acts_by_obj[act]
                for p in model.par_dmr_subscrs_by_act[act_uindx]:
                    model.c1c.add(model.var_act_indic[act_uindx] >= model.var_dmr_indic[p]) 

        def c2_rule( model,p):
            return model.par_dmr_dv[p]*model.var_dmr_utilization[p] >= model.par_min_as_route_dv*model.var_dmr_indic[p]
        model.c2 =pe.Constraint ( model.dmrs,  rule=c2_rule)

        def c3_rule( model,a):
            return model.var_act_indic[a] >=  model.var_activity_utilization[a]
        model.c3 =pe.Constraint ( model.acts,  rule=c3_rule)  

        def c3c_rule( model,p):
            return model.var_dmr_indic[p] >=  model.var_dmr_utilization[p]
        model.c3c =pe.Constraint ( model.dmrs,  rule=c3c_rule)

        # keep track of this, so we know to warn user about default transition time usage
        used_default_transition_time = False

        #  intra-satellite activity overlap constraints [4],[5],[5b]
        model.c4  = pe.ConstraintList()
        model.c5  = pe.ConstraintList()
        model.c5b  = pe.ConstraintList()
        for sat_indx in range (self.num_sats):
            num_sat_acts = len(sat_acts[sat_indx])
            for  first_act_indx in  range (num_sat_acts):

                # from circinus_tools import debug_tools
                # debug_tools.debug_breakpt()

                act1 = sat_acts[sat_indx][first_act_indx]
                # get the unique index into model.acts
                act1_uindx = all_acts_by_obj[act1]
                length_1 = model.par_act_dv[act1_uindx]*model.var_activity_utilization[act1_uindx]/act1.ave_data_rate
                model.c5b.add( length_1 >= model.var_act_indic[act1_uindx] * self.min_act_duration_s[type(act1)])
                
                for  second_act_indx in  range (first_act_indx+1,num_sat_acts):
                    act2 = sat_acts[sat_indx][ second_act_indx]
                    # get the unique index into model.acts
                    act2_uindx = all_acts_by_obj[act2]

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

                    #  if the activities overlap in center time, then it's not possible to have sufficient transition time between them
                    #  add constraint to rule out the possibility of scheduling both of them
                    elif (act2.center - act1.center).total_seconds() <= transition_time_req:
                        model.c4.add( model.var_act_indic[act1_uindx]+ model.var_act_indic[act2_uindx] <= 1)

                    # If they don't overlap in center time, but they do overlap to some amount, then we need to constrain their end and start times to be consistent with one another
                    else:
                        center_time_diff = (act2.center - act1.center).total_seconds()
                        # this is the adjustment added to the center time to get to the start or end of the activity
                        time_adjust_1 = model.par_act_dv[act1_uindx]*model.var_activity_utilization[act1_uindx]/2/act1.ave_data_rate
                        time_adjust_2 = model.par_act_dv[act2_uindx]*model.var_activity_utilization[act2_uindx]/2/act2.ave_data_rate
                        constr_disable_1 = self.big_M_act_t_dur_s*(1-model.var_act_indic[act1_uindx])
                        constr_disable_2 = self.big_M_act_t_dur_s*(1-model.var_act_indic[act2_uindx])
                        model.c5.add( center_time_diff - time_adjust_1 - time_adjust_2 + constr_disable_1 + constr_disable_2 >= transition_time_req)


        # inter-satellite downlink overlap constraints [9],[10]
        model.c9  = pe.ConstraintList()
        model.c10  = pe.ConstraintList()
        for sat_indx in range (self.num_sats):
            num_sat_acts = len(sat_dlnks[sat_indx])
            
            for other_sat_indx in range (self.num_sats):
                if other_sat_indx == sat_indx:
                    continue

                num_other_sat_acts = len(sat_dlnks[other_sat_indx])

                for  sat_act_indx in  range (num_sat_acts):

                    # from circinus_tools import debug_tools
                    # debug_tools.debug_breakpt()

                    act1 = sat_dlnks[sat_indx][sat_act_indx]
                    # get the unique index into model.acts
                    act1_uindx = all_acts_by_obj[act1]
                    
                    for  other_sat_act_indx in  range (num_other_sat_acts):
                        act2 = sat_dlnks[other_sat_indx][other_sat_act_indx]
                        # get the unique index into model.acts
                        act2_uindx = all_acts_by_obj[act2]

                        assert(type(act1) == DlnkWindow and type(act2) == DlnkWindow)

                        # todo: include!
                        # # if they're not looking at the same gs, then constraints don't apply
                        # if not act1.gs_indx == act2.gs_indx:
                        #     continue

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

                        #  if the activities overlap in center time, then it's not possible to have sufficient transition time between them
                        #  add constraint to rule out the possibility of scheduling both of them
                        elif (act2.center - act1.center).total_seconds() <= transition_time_req:
                            model.c9.add( model.var_act_indic[act1_uindx]+ model.var_act_indic[act2_uindx] <= 1)

                        # If they don't overlap in center time, but they do overlap to some amount, then we need to constrain their end and start times to be consistent with one another
                        else:
                            center_time_diff = (act2.center - act1.center).total_seconds()
                            # this is the adjustment added to the center time to get to the start or end of the activity
                            time_adjust_1 = model.par_act_dv[act1_uindx]*model.var_activity_utilization[act1_uindx]/2/act1.ave_data_rate
                            time_adjust_2 = model.par_act_dv[act2_uindx]*model.var_activity_utilization[act2_uindx]/2/act2.ave_data_rate
                            constr_disable_1 = self.big_M_act_t_dur_s*(1-model.var_act_indic[act1_uindx])
                            constr_disable_2 = self.big_M_act_t_dur_s*(1-model.var_act_indic[act2_uindx])
                            model.c10.add( center_time_diff - time_adjust_1 - time_adjust_2 + constr_disable_1 + constr_disable_2 >= transition_time_req)

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()


        if verbose:
            if used_default_transition_time:
                print('\nWarning: used default transition time for inter- or intra- satellite activity timing constraints\n')


        #  energy constraints [6]
        model.c6  = pe.ConstraintList()
        for sat_indx in range (self.num_sats): 

            # tp_indx serves as an index into the satellite activity dance cards
            for tp_indx in model.es_timepoint_indcs:
                #  constraining first time step to initial energy storage
                #  continue for loop afterwards because no geq/leq constraints needed for this index
                if tp_indx == 0:
                    model.c6.add( model.var_sats_estore[sat_indx,tp_indx] ==  model.par_sats_estore_initial[sat_indx])
                    continue 

                #  minimum and maximum storage constraints
                model.c6.add( model.var_sats_estore[sat_indx,tp_indx] >= model.par_sats_estore_min[sat_indx])
                model.c6.add( model.var_sats_estore[sat_indx,tp_indx] <= model.par_sats_estore_max[sat_indx])

                if self.enforce_energy_storage_constr:
                    # determine activity energy consumption
                    charging = True
                    activity_delta_e = 0 
                    #  get the activities that were active during the time step immediately preceding time point
                    activities = es_act_dancecards[sat_indx].get_objects_at_ts_pre_tp_indx(tp_indx)
                    # activities may be none if nothing is happening at timestep, to minimize RAM usage
                    if activities:
                        for act in activities:
                            #  if this is a "standard activity" that we can choose to perform or not
                            if type(act) in self.standard_activities:
                                act_uindx = all_acts_by_obj[act]
                                activity_delta_e += (
                                    model.par_sats_edot_by_act[sat_indx][self.act_type_map[type(act)]] 
                                    * model.var_activity_utilization[act_uindx]
                                    * model.par_resource_delta_t
                                )

                            if type(act) == XlnkWindow:
                                act_uindx = all_acts_by_obj[act]

                                if self.use_symmetric_xlnk_windows:
                                    xlnk_edot = model.par_sats_edot_by_act[sat_indx]['xlnk-tx']
                                # note that in the case where we're not using symmetric xlnk windows, both satellites have a copy of each unidirectional window in their list of activity windows, so it'll be added appropriately at each unique sat_indx
                                else:
                                    is_rx = act.is_rx(sat_indx)
                                    if is_rx: 
                                        xlnk_edot = model.par_sats_edot_by_act[sat_indx]['xlnk-rx']
                                    else: 
                                        xlnk_edot = model.par_sats_edot_by_act[sat_indx]['xlnk-tx']

                                activity_delta_e += (
                                    xlnk_edot 
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


        #  data storage constraints [7]
        model.c7  = pe.ConstraintList()
        for sat_indx in range (self.num_sats): 

            # tp_indx serves as an index into the satellite data storage dance cards
            for tp_indx in model.ds_timepoint_indcs:

                # todo: add in an intital data volume value?
                #  maximum storage constraints
                model.c7.add( model.var_sats_dstore[sat_indx,tp_indx] <= model.par_sats_dstore_max[sat_indx])

                if self.enforce_data_storage_constr:
                    routes_storing = ds_route_dancecards[sat_indx][tp_indx]

                    # todo: this may be too slow to use == below. Change it back to >= and use a better approach to extract real data storage values in extract_resource_usage below?
                    # constrain the minimum data storage at this tp by the amount of data volume being buffered by each route passing through the sat
                    if routes_storing:
                        model.c7.add( model.var_sats_dstore[sat_indx,tp_indx] == sum(model.par_dmr_dv[p]*model.var_dmr_utilization[p] for p in routes_storing))
                    else:
                        model.c7.add( model.var_sats_dstore[sat_indx,tp_indx] == 0)



        #  observation latency score factor constraints [8]
        model.c8  = pe.ConstraintList()
        for obs_indx in model.obs_acts:
            dmrs_obs = model.par_dmr_subscrs_by_obs_act[obs_indx]

            #  sort the latency score factors for all the dmrs for this observation in increasing order -  important for constraint construction
            dmrs_obs.sort(key= lambda p: dmr_latency_sf_by_dmr_indx[p])

            num_dmrs_obs = len(dmrs_obs)
            #  initial constraint -  score factor for this observation will be equal to zero if no dmrs for this obs were chosen
            model.c8.add( model.var_dmr_latency_sf_by_obs_indx[obs_indx] <= 0 + self.big_M_lat * sum(model.var_dmr_indic[p] for p in dmrs_obs) )
            
            for dmr_obs_indx in range(num_dmrs_obs):
                #  add constraint that score factor for observation is less than or equal to the score factor for this dmr_obs_indx, plus any big M terms for any dmrs with larger score factors.
                #  what this does is effectively disable the constraint for the score factor for this dmr_obs_indx if any higher score factor dmrs were chosen 
                model.c8.add( model.var_dmr_latency_sf_by_obs_indx[obs_indx] <= 
                    dmr_latency_sf_by_dmr_indx[dmrs_obs[dmr_obs_indx]] + 
                    self.big_M_lat * sum(model.var_dmr_indic[p] for p in dmrs_obs[dmr_obs_indx+1:num_dmrs_obs]) )

                #  note: use model.c8[indx].expr.to_string()  to print out the constraint in a human readable form
                #                ^ USES BASE 1 INDEXING!!! WTF??
                

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        ##############################
        #  Make objective
        ##############################


        #  determine which time points to use as "spot checks" on resource margin. These are the points that will be used in the objective function for maximizing resource margin
        timepoint_spacing = ceil(es_num_timepoints/self.resource_margin_obj_num_timepoints)
        # need to turn the generator into a list for slicing
        #  note: have to get the generator again
        decimated_tp_indcs = list(self.es_time_getter_dc.get_tp_indcs())[::timepoint_spacing]
        rsrc_norm_f = len(decimated_tp_indcs) * len(model.sats)

        def obj_rule(model):
            return (
                self.obj_weights['obs_dv'] * 1/model.par_total_obs_dv * sum(model.par_dmr_dv[p]*model.var_dmr_utilization[p] for p in model.dmrs) 
                # + self.obj_weights['route_latency'] * 1/len(model.obs_acts) * sum(dmr_latency_sf_by_dmr_indx[p]*model.var_dmr_utilization[p] for p in model.dmrs)
                + self.obj_weights['route_latency'] * 1/len(model.obs_acts) * sum(model.var_dmr_latency_sf_by_obs_indx[o] for o in model.obs_acts)
                + self.obj_weights['energy_storage'] * 1/rsrc_norm_f * sum(model.var_sats_estore[sat_indx,tp_indx]/model.par_sats_estore_max[sat_indx] for tp_indx in decimated_tp_indcs for sat_indx in model.sats)
            )
        model.obj = pe.Objective( rule=obj_rule, sense=pe.maximize )

        self.model = model

    # taken in part from Jeff Menezes' code at https://github.mit.edu/jmenezes/Satellite-MILP/blob/master/sat_milp_pyomo.py
    def solve(self):

        solver = po.SolverFactory(self.solver_name)
        if self.solver_name == 'gurobi':
            # note default for this is 1e-4, or 0.01%
            solver.options['TimeLimit'] = self.solver_params['max_runtime_s']
            solver.options['MIPGap'] = self.solver_params['optimality_gap']
            solver.options['IntFeasTol'] = self.solver_params['integer_feasibility_tolerance']
        elif self.solver_name == 'cplex':
            solver.options['timelimit'] = self.solver_params['max_runtime_s']
            solver.options['mip_tolerances_mipgap'] = self.solver_params['optimality_gap']
            solver.options['mip_tolerances_integrality'] = self.solver_params['integer_feasibility_tolerance']


            solver.options['TimeLimit'] = 100 #self.solver_params['gurobi']['max_runtime_s']
            solver.options['MIPGap'] = self.solver_params['gurobi']['optimality_gap']
            solver.options['IntFeasTol'] = self.solver_params['gurobi']['integer_feasibility_tolerance']
            # solver.options['Cuts'] = 0
            # solver.options['MIPFocus'] = 1 #for finding feasible solutions quickly
            # solver.options['MIPFocus'] = 3 #for lowering the mip gap

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

    def print_sol(self):
        for v in self.model.component_objects(pe.Var, active=True):
            if str (v) =='var_activity_utilization': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" %d %0.3f  %s"%(index, val, self.all_acts_by_indx[index]))

            elif str (v) =='var_dmr_utilization': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" ",index, val)
            
            elif str (v) =='var_dmr_indic': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" ",index, val)

    def extract_utilized_routes( self, copy_routes = True, verbose = False):
        # scheduled_dv
        scheduled_routes_flat = []

        if verbose:
            print ('utilized routes:')

        # commenting out for now because I'm not sure this is actually useful, considering we already have to modify the windows below...
        # def copy_choice(route):
        #     if copy_routes:
        #         # this will copy everything but the underlying activity windows. 
        #         return copy(route)
        #     else:
        #         return route

        # figure out which dmrs were used, and add the scheduled data volume to each
        for p in self.model.dmrs:
            if pe.value(self.model.var_dmr_indic[p]) >= 1.0 - self.binary_epsilon:
                # scheduled_route =  copy_choice (self.routes_flat[p]) 
                scheduled_route =  self.routes_flat[p]
                scheduled_route.set_scheduled_dv_frac(pe.value(self.model.var_dmr_utilization[p]))  #* self.model.par_dmr_dv[p])
                scheduled_routes_flat. append (scheduled_route)

                if verbose:
                    print(scheduled_route)

        # examine the schedulable data volume for every activity window, checking as we go that the data volume is sufficient for at least the route in which the window is found
        #  note that this code is slightly inefficient because it might duplicate windows across routes. that's fine though, because we're thorough in checking across all routes
        # note: dmr is for DataMultiRoute
        wind_sched_dv_check = {}
        for dmr in scheduled_routes_flat:
            # wind may get set multiple times due to Windows appearing across routes, but that's not really a big deal
            for wind in dmr.get_winds():
                act_indx = self.all_acts_by_obj[wind]
                # wind_sched_dv_check[wind] = wind.data_vol * pe.value(self.model.var_act_indic[act_indx])
                wind_sched_dv_check[wind] = wind.data_vol * pe.value(self.model.var_activity_utilization[act_indx])
                #  initialize this while we're here
                wind.scheduled_data_vol = 0

            # similar to winds, set dr data vols to 0
            for dr in dmr.data_routes:
                # dmrs should not be sharing drs across themselves, so no scheduled dv should be seen yet
                assert(dr.scheduled_dv == const.UNASSIGNED)
                dr.scheduled_dv = dmr.scheduled_dv_by_dr[dr] 


        #  now we want to mark the real scheduled data volume for every window. We need to do this separately because the model.var_act_indic continuous variables only give an upper bound on the data volume for an activity. we only actually need to use as much data volume as the data routes want to push through the window
        #  add data volume for every route passing through every window
        for dmr in scheduled_routes_flat:
            for wind in dmr.get_winds():
                wind.scheduled_data_vol += dmr.scheduled_dv_for_wind(wind)


        # update the window beginning and end times based upon their amount of scheduled data volume
        # keep track of which ones we've updated, because we should only update once
        updated_winds = set()
        for dmr in scheduled_routes_flat:
            # validate the data multi route (and in turn, the scheduled data vols of all the data routes under it)
            dmr.validate()

            for wind in dmr.get_winds():
                #  this check should be at least as big as the scheduled data volume as calculated from all of the route data volumes. (it's not constrained from above, so it could be bigger)
                if wind_sched_dv_check[wind] < wind.scheduled_data_vol - self.dv_epsilon:
                    raise RuntimeWarning('inconsistent activity scheduling results: activity window data volume determined from route dvs (%f) is greater than dv from var_act_indic (%f) [dmr: %s, wind: %s]'%(wind.scheduled_data_vol,wind_sched_dv_check[wind],dmr,wind))

                if not wind in updated_winds:
                    # note that the line below seems like it may break the scheduled times for activities by specifying a minimum activity duration. however, this minimum activity duration is already accounted for in scheduling constraints
                    wind.update_duration_from_scheduled_dv (min_duration_s=self.min_act_duration_s[type(wind)])
                    updated_winds.add(wind)

        return scheduled_routes_flat

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


