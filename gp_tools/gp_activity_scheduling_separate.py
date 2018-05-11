# contains model and Solver for global planner activity scheduling capability
# 
# @author Kit Kennedy
#

from  datetime import timedelta
from copy import  deepcopy
from math import ceil
from collections import OrderedDict

from pyomo import environ  as pe
from pyomo import opt  as po

import numpy as np

from circinus_tools  import time_tools as tt
from .gp_activity_scheduling_super import  GPActivityScheduling
from circinus_tools  import  constants as const
from circinus_tools.scheduling.custom_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow
from circinus_tools.scheduling.schedule_objects import Dancecard
from circinus_tools.scheduling.routing_objects import DataMultiRoute

def print_verbose(string,verbose=False):
    if verbose:
        print(string)

class GPActivitySchedulingSeparate(GPActivityScheduling):
    """docstring for GP activity scheduling"""
    
    def __init__(self,gp_params):
        """initializes based on parameters
        
        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """

        super().__init__(gp_params)

        # number of routes found after filtering
        self.num_routes_filt = const.UNASSIGNED

    def get_stats(self,verbose=True):

        num_winds_per_route = [len(dmr.get_winds()) for dmr in self.routes_filt]
        num_routes_by_act = {act:len(self.dmr_ids_by_act_windid[act_indx]) for act,act_indx in self.all_act_windids_by_obj.items()}

        stats = {}
        stats['num_routes'] = sum([len ( self.routes_filt)])
        stats['num_acts'] = sum([len ( self.all_acts_windids)])
        stats['num_obs_acts'] = sum([len ( self.obs_windids)])
        stats['num_link_acts'] = sum([len ( self.lnk_windids)])
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
        num_routes_xlnk = [num for act,num in num_routes_by_act.items() if type(act) == XlnkWindow]
        stats['ave_num_routes_per_xlnk'] = np.mean(num_routes_xlnk) if len(num_routes_xlnk) > 0 else const.UNASSIGNED
        stats['max_num_routes_per_xlnk'] = np.max(num_routes_xlnk) if len(num_routes_xlnk) > 0 else const.UNASSIGNED

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

    def get_act_model_objs(self,act):

        model_objs_act = {
            'act_object': act,
            'var_dv_utilization': self.model.var_activity_utilization[act.window_ID]*self.model.par_act_capacity[act.window_ID],
            'par_dv_capacity': self.model.par_act_capacity[act.window_ID],
            'var_act_indic': self.model.var_act_indic[act.window_ID],
        }

        return model_objs_act

    def get_activity_structs( self,routes_filt):

        #  all activities are uniquely indexed. these structures keep track of those, and the mapping to activity objects
        all_acts_windids = []
        #  these structures are for lookup in both directions
        all_acts_by_windid = {}
        all_act_windids_by_obj = {}

        # these structures keep track of the subset of unique indices that correspond to observations and links. Also we keep track of what data routes correspond to an activity
        all_lnk_windids = []
        all_obs_windids = []
        # dmr is data multi-route
        dmr_ids_by_link_act_windid = {}
        dmr_ids_by_obs_act_windid = {}
        dmr_ids_by_act_windid = {}
        dv_by_link_act_windid = {}
        dv_by_obs_act_windid = {}
        dv_by_act_windid = {}

        sats_acts = [[] for sat_indx in range (self.num_sats)]
        sats_dlnks = [[] for sat_indx in range (self.num_sats)]

        for dmr in routes_filt:
            for act in dmr.get_winds():

                # We have already filtered routes, but we also need to filter activities because there may be activities from existing routes that are outside of the planning window. we cannot fully enforce constraints on such windows, so we need to drop them
                if act.start < self.planning_start_dt or act.end > self.planning_end_dt:
                    continue

                # if we haven't yet seen this activity, then add it to bookkeeping
                if not act in all_act_windids_by_obj.keys():
                    # use activity window ID as a unique id for window in the model
                    act_windid = act.window_ID
                    all_acts_windids.append(act_windid)
                    all_acts_by_windid[act_windid] = act
                    all_act_windids_by_obj[act] = act_windid
                    dmr_ids_by_act_windid[act_windid] = []
                    dmr_ids_by_act_windid[act_windid].append (dmr.ID)
                    dv_by_act_windid[act_windid] = act.data_vol

                    # also need to add it to the list and dictionary for observations
                    if type(act) == ObsWindow:
                        sats_acts[act.sat_indx].append(act)
                        all_obs_windids.append(act_windid)
                        dmr_ids_by_obs_act_windid[act_windid] = []
                        dmr_ids_by_obs_act_windid[act_windid].append (dmr.ID)
                        dv_by_obs_act_windid[act_windid] = act.data_vol

                    # also need to add it to the list and dictionary for links
                    if type(act) == DlnkWindow:
                        all_lnk_windids.append(act_windid)
                        dmr_ids_by_link_act_windid[act_windid] = []
                        dmr_ids_by_link_act_windid[act_windid].append (dmr.ID)
                        dv_by_link_act_windid[act_windid] = act.data_vol
                        sats_acts[act.sat_indx].append(act)
                        # grab the dlnks for each sat too, while we're looping through
                        sats_dlnks[act.sat_indx].append(act)

                    if type(act) == XlnkWindow:
                        all_lnk_windids.append(act_windid)
                        dmr_ids_by_link_act_windid[act_windid] = []
                        dmr_ids_by_link_act_windid[act_windid].append (dmr.ID)
                        dv_by_link_act_windid[act_windid] = act.data_vol
                        sats_acts[act.sat_indx].append(act)
                        sats_acts[act.xsat_indx].append(act)

                #  if we have already seen the activity,  then just need to update the appropriate structures
                else:
                    act_windid = all_act_windids_by_obj[act]
                    dmr_ids_by_act_windid[act_windid].append (dmr.ID)

                    # add the data route index
                    if type(act) == ObsWindow:
                        dmr_ids_by_obs_act_windid[act_windid].append (dmr.ID)
                    if type(act) == DlnkWindow or type(act) == XlnkWindow:
                        dmr_ids_by_link_act_windid[act_windid].append (dmr.ID)

        # print (dmr_ids_by_link_act_windid)
        # print (dmr_ids_by_obs_act_windid)

        #  sort the activities, because we'll need that for constructing constraints
        for sat_indx in range (self.num_sats):
            sats_acts[sat_indx].sort(key=lambda x: x.center)
            sats_dlnks[sat_indx].sort(key=lambda x: x.center)

        return sats_acts,sats_dlnks,all_acts_windids,dmr_ids_by_act_windid,dv_by_act_windid,all_act_windids_by_obj,all_acts_by_windid,all_obs_windids,dmr_ids_by_obs_act_windid,dv_by_obs_act_windid,all_lnk_windids,dmr_ids_by_link_act_windid,dv_by_link_act_windid
                    

    def filter_routes( self,new_routes,existing_routes):

        routes_filt = []

        def process_dmr(dmr,start_filt_dt,filter_opt):
            dmr_start = dmr.get_obs().start
            dmr_end = dmr.get_dlnk().end

            # check if all act windows in route are completely within the planning window. Pass if not.
            if filter_opt=='totally_within' and (dmr_start < start_filt_dt or dmr_end > self.planning_end_dt):
                pass
            # check if at least one act window in route is partially within the planning window. Pass if not.
            elif filter_opt=='partially_within' and (dmr_end < start_filt_dt or dmr_start > self.planning_end_dt):
                pass
            elif dmr.get_dlnk().duration.total_seconds() < self.min_act_duration_s[DlnkWindow]:
                print('discarding too short dlnk window')
                pass
            else:
                routes_filt.append (dmr)
                return True

            return False
            
        # want new routes to be entirely within non-fixed planning window (otherwise we could miss enforcing all of the constraints imposed by the route)
        for dmr in new_routes:
            process_dmr(dmr,self.planning_fixed_end_dt,"totally_within")

        # want existing routes to be at least partially within whole planning window, so that we will be able to later filter for all act windows on which they impose constraints (but note that we don't have to worry anymore about the constraints on the route, because the input utilization number for the route accounts for it) 
        for dmr in existing_routes:
            added = process_dmr(dmr,self.planning_start_dt,"partially_within")

        return routes_filt

    def get_dmr_latency_score_factors(self,routes_by_dmr_id,dmr_ids_by_obs_act_windid):

        latency_sf_by_dmr_id = {}

        # loop through all satellites and their downlinks
        for obs,dmr_ids in dmr_ids_by_obs_act_windid.items():
            latencies = []
            for dmr_id in dmr_ids:
                dmr = routes_by_dmr_id[dmr_id]

                latencies.append(
                    dmr.get_latency(
                        'minutes',
                        obs_option = self.latency_params['obs'], 
                        dlnk_option = self.latency_params['dlnk']
                    )
                )

             #  the shortest latency dmr (DataMultiRoute) for this observation has a score factor of 1.0, and the score factors for the other dmrs decrease as the inverse of increasing latency
            min_lat = min(latencies)
            min_lat = max(min_lat, self.min_latency_for_sf_1_mins)
            for lat_indx, lat in enumerate(latencies):
                dmr_id = dmr_ids[lat_indx]
                lat = max(lat, self.min_latency_for_sf_1_mins)
                latency_sf_by_dmr_id[dmr_id] = min_lat/lat
 
        return latency_sf_by_dmr_id

    def make_model ( self,new_routes, existing_routes, utilization_by_existing_route_id, ecl_winds, verbose = True):
        # important assumption: all activity window IDs are unique!

        model = pe.ConcreteModel()
        self.model = model

        # filter the routes to make sure that  none of their activities fall outside the scheduling window
        routes_filt = self.filter_routes(new_routes,existing_routes)
        routes_by_dmr_id = {dmr.ID:dmr for dmr in routes_filt}

        self.routes_filt = routes_filt
        fixed_routes_ids = {dmr.ID for dmr in existing_routes}
        self.utilization_by_existing_route_id = utilization_by_existing_route_id
        self.routes_by_dmr_id = routes_by_dmr_id
        self.ecl_winds = ecl_winds

        self.num_routes_filt = len(routes_filt)
        if self.num_routes_filt == 0:
            if verbose:
                print('No routes found! Quitting separate AS early')
            self.model_constructed = False
            return

        #  should only be using data multi-route objects for activity scheduling, even if they're just a shallow wrapper around a DataRoute
        for dmr in routes_filt:
            assert(type(dmr) == DataMultiRoute)

        ##############################
        #  Make indices/ subscripts
        ##############################

        try:
            (sats_acts,
                sats_dlnks,
                all_acts_windids,
                dmr_ids_by_act_windid,
                dv_by_act_windid,
                all_act_windids_by_obj,
                all_acts_by_windid,
                all_obs_windids,
                dmr_ids_by_obs_act_windid,
                dv_by_obs_act_windid,
                all_lnk_windids,
                dmr_ids_by_link_act_windid,
                dv_by_link_act_windid) =  self.get_activity_structs(routes_filt)

            latency_sf_by_dmr_id =  self.get_dmr_latency_score_factors(
                routes_by_dmr_id,
                dmr_ids_by_obs_act_windid
            )

            # verify that all acts found are within the planning window, otherwise we may end up with strange results
            for sat_acts in sats_acts:
                for act in sat_acts:
                    if act.start < self.planning_start_dt or act.end > self.planning_end_dt:
                        raise RuntimeWarning('Activity is out of planning window range (start %s, end %s): %s'%(self.planning_start_dt,self.planning_end_dt,act))

            # construct a set of dance cards for every satellite, 
            # each of which keeps track of all of the activities of satellite 
            # can possibly execute at any given time slice delta T. 
            # this is for constructing energy storage constraints
            # using resource_delta_t_s because this dancecard is solely for use in constructing resource constraints
            # note that these dancecards will baloon in size pretty quickly as the planning_end_dt increases. However most of the complexity they introduce is before planning_end_obs_xlnk_dt because that's the horizon where obs,xlnk actitivities are included. After that there should only be sparse downlinks
            es_act_dancecards = [Dancecard(self.planning_start_dt,self.planning_end_dt,self.resource_delta_t_s,item_init=None,mode='timestep') for sat_indx in range (self.num_sats)]
            
            for sat_indx in range (self.num_sats): 
                es_act_dancecards[sat_indx].add_winds_to_dancecard(sats_acts[sat_indx])
                es_act_dancecards[sat_indx].add_winds_to_dancecard(ecl_winds[sat_indx])

            # this is for data storage
            # for each sat/timepoint, we store a list of all those data multi routes that are storing data on the sat at that timepoint
            ds_route_dancecards = [Dancecard(self.planning_start_dt,self.planning_end_dt,self.resource_delta_t_s,item_init=None,mode='timepoint') for sat_indx in range (self.num_sats)]
            
            # add data routes to the dancecard
            for dmr in routes_filt:
                # list of type routing_objects.SatStorageInterval
                dmr_ds_intervs = dmr.get_data_storage_intervals()

                for interv in dmr_ds_intervs:
                    # store the dmr object at this timepoint
                    ds_route_dancecards[interv.sat_indx].add_item_in_interval(dmr.ID,interv.start,interv.end)

        except IndexError:
            raise RuntimeWarning('sat_indx out of range. Are you sure all of your input files are consistent? (including pickles)')        
        self.dmr_ids_by_act_windid = dmr_ids_by_act_windid
        self.all_acts_windids = all_acts_windids
        self.obs_windids = all_obs_windids
        self.lnk_windids = all_lnk_windids
        self.all_act_windids_by_obj = all_act_windids_by_obj
        self.all_acts_by_windid = all_acts_by_windid

        # these dmr subscripts probably should've been done using the unique IDs for the objects, rather than their arbitrary locations within a list. Live and learn, hélas...

        #  subscript for each dmr (data multi route) p  (p index is a hold-over from when I called them paths)
        model.dmr_ids = pe.Set(initialize= routes_by_dmr_id.keys())
        #  subscript for each activity a
        model.act_windids = pe.Set(initialize= all_acts_windids)
        #  subscript for each satellite
        model.sat_indcs = pe.Set(initialize=  range ( self.num_sats))

        # timepoints is the indices, which starts at 0 
        #  NOTE: we assume the same time system for every satellite
        self.es_time_getter_dc = es_act_dancecards[0]
        es_num_timepoints = es_act_dancecards[0].num_timepoints
        model.es_timepoint_indcs = pe.Set(initialize=  self.es_time_getter_dc.get_tp_indcs())

        self.ds_time_getter_dc = ds_route_dancecards[0]
        model.ds_timepoint_indcs = pe.Set(initialize=  self.ds_time_getter_dc.get_tp_indcs())

        #  unique indices for observation and link acts
        model.obs_windids = pe.Set(initialize= all_obs_windids)
        model.lnk_windids = pe.Set(initialize= all_lnk_windids)

        if self.solver_name == 'gurobi' or self.solver_name == 'cplex':
            int_feas_tol = self.solver_params['integer_feasibility_tolerance']
        else:
            raise NotImplementedError

        for p,lat_sf in latency_sf_by_dmr_id.items():        
            if lat_sf > int_feas_tol*self.big_M_lat:
                raise RuntimeWarning('big_M_lat (%f) is not large enough for latency score factor %f and integer feasibility tolerance %f (dmr index %d)'%(self.big_M_lat,lat_sf,int_feas_tol,p))

        # for act_obj in all_act_windids_by_obj.keys():
        #     if 2*(act_obj.end-act_obj.start).total_seconds() > self.big_M_act_t_dur_s:
        #         raise RuntimeWarning('big_M_act_t_dur_s (%f) is not large enough for act of duration %s and integer feasibility tolerance %f (act string %s)'%(self.big_M_act_t_dur_s,act_obj.end-act_obj.start,int_feas_tol,act_obj))
        #     # if 2*(act_obj.end-act_obj.start).total_seconds() > (1-int_feas_tol) *  self.big_M_act_t_dur_s:
        #         # raise RuntimeWarning('big_M_act_t_dur_s (%f) is not large enough for act of duration %s and integer feasibility tolerance %f (act string %s)'%(self.big_M_act_t_dur_s,act_obj.end-act_obj.start,int_feas_tol,act_obj))

        #     if 2*act_obj.data_vol > (1-int_feas_tol) * self.big_M_act_dv:
        #         raise RuntimeWarning('big_M_act_dv (%f) is not large enough for act of dv %f and integer feasibility tolerance %f (act string %s)'%(self.big_M_act_dv,act_obj.data_vol,int_feas_tol,act_obj))
                


        ##############################
        #  Make parameters
        ##############################

        model.par_min_obs_dv_dlnk_req = pe.Param (initialize=self.min_obs_dv_dlnk_req)
        model.par_obs_capacity = pe.Param(model.obs_windids,initialize =dv_by_obs_act_windid)
        model.par_total_obs_dv = sum(dv_by_obs_act_windid.values())
        model.par_link_capacity = pe.Param(model.lnk_windids,initialize =dv_by_link_act_windid)
        model.par_act_capacity = pe.Param(model.act_windids,initialize =dv_by_act_windid)
        #  data volume for each data multi-route
        model.par_dmr_dv = pe.Param(model.dmr_ids,initialize ={ dmr.ID: dmr.data_vol for dmr in routes_filt})
        #  data volume for each activity in each data multi-route

        # todo: will probably have to fix this to not be the full set multiplication of model.dmr_ids, model.act_windids - can use dmr_ids_by_act_windid
        model.par_dmr_act_dv = pe.Param(
            model.dmr_ids,
            model.act_windids,
            initialize = { (dmr.ID,self.all_act_windids_by_obj[act]): 
                dmr.data_vol_for_wind(act) for dmr in routes_filt for act in dmr.get_winds()
            }
        )

        # each of these is essentially a dictionary indexed by link or obs act indx, with  the value being a list of dmr indices that are included within that act
        # these are valid indices into model.dmr_ids
        model.par_dmr_subscrs_by_link_act = pe.Param(model.lnk_windids,initialize =dmr_ids_by_link_act_windid)
        model.par_dmr_subscrs_by_obs_act = pe.Param(model.obs_windids,initialize =dmr_ids_by_obs_act_windid)
        model.par_dmr_subscrs_by_act = pe.Param(model.act_windids,initialize =dmr_ids_by_act_windid)


        if self.energy_unit == "Wh":
            model.par_resource_delta_t = pe.Param (initialize= self.resource_delta_t_s/3600)
        else: 
            raise NotImplementedError
        model.par_sats_estore_initial = pe.Param ( model.sat_indcs,initialize= { i: item for i,item in enumerate (self.sats_init_estate_Wh)})
        model.par_sats_estore_min = pe.Param ( model.sat_indcs,initialize= { i: item for i,item in enumerate (self.sats_emin_Wh)})
        model.par_sats_estore_max = pe.Param ( model.sat_indcs,initialize= { i: item for i,item in enumerate (self.sats_emax_Wh)})
        model.par_sats_edot_by_mode = pe.Param ( model.sat_indcs,initialize= { i: item for i,item in enumerate (self.sats_edot_by_mode_W)})

        model.par_sats_dstore_min = pe.Param ( model.sat_indcs,initialize= { i: item for i,item in enumerate (self.sats_dmin_Mb)})
        model.par_sats_dstore_max = pe.Param ( model.sat_indcs,initialize= { i: item for i,item in enumerate (self.sats_dmax_Mb)})

        ##############################
        #  Make variables
        ##############################


        # activity utilization variable indicating how much of an activity's capacity is used [1]
        model.var_activity_utilization  = pe.Var (model.act_windids, bounds =(0,1))
        # dmr utilization variable indicating how much of a dmr's capacity is used [2]
        model.var_dmr_utilization  = pe.Var (model.dmr_ids, bounds =(0,1))
        #  indicator variables for whether or not dmrs [3] and activities [4] have been chosen
        model.var_dmr_indic  = pe.Var (model.dmr_ids, within = pe.Binary)
        # model.var_dmr_indic  = pe.Var (model.dmr_ids, bounds =(0,1))
        model.var_act_indic  = pe.Var (model.act_windids, within = pe.Binary)
        # model.var_act_indic  = pe.Var (model.act_windids, bounds =(0,1))

        
        # satellite energy storage
        model.var_sats_estore  = pe.Var (model.sat_indcs,  model.es_timepoint_indcs,  within = pe.NonNegativeReals)

        # satellite data storage (data buffers)
        model.var_sats_dstore  = pe.Var (model.sat_indcs,  model.ds_timepoint_indcs,  within = pe.NonNegativeReals)

        model.var_latency_sf_obs = pe.Var (model.obs_windids,  bounds = (0,1.0))
        
        allow_act_timing_constr_violations = False
        if allow_act_timing_constr_violations:
            print('allow_act_timing_constr_violations is True')

        #  variables for handling the allowance of inter-activity timing constraint violations. these are only generated if allow_act_timing_constr_violations is True
        model.var_intra_sat_act_constr_violations = pe.VarList()
        model.var_inter_sat_act_constr_violations = pe.VarList()
        model.intra_sat_act_constr_bounds  = pe.ConstraintList()
        model.inter_sat_act_constr_bounds  = pe.ConstraintList()

        #  stores all of the lower bounds of the constraint violation variables, for use in normalization for objective function
        min_var_intra_sat_act_constr_violation_list = [] 
        min_var_inter_sat_act_constr_violation_list = [] 

        constraint_violation_model_objs = {}
        constraint_violation_model_objs['intra_sat_act_constr_violation_acts_list'] = []
        constraint_violation_model_objs['inter_sat_act_constr_violation_acts_list'] = []
        constraint_violation_model_objs['var_intra_sat_act_constr_violations'] = model.var_intra_sat_act_constr_violations
        constraint_violation_model_objs['var_inter_sat_act_constr_violations'] = model.var_inter_sat_act_constr_violations
        constraint_violation_model_objs['intra_sat_act_constr_bounds'] = model.intra_sat_act_constr_bounds
        constraint_violation_model_objs['inter_sat_act_constr_bounds'] = model.inter_sat_act_constr_bounds
        constraint_violation_model_objs['min_var_intra_sat_act_constr_violation_list'] = min_var_intra_sat_act_constr_violation_list 
        constraint_violation_model_objs['min_var_inter_sat_act_constr_violation_list'] = min_var_inter_sat_act_constr_violation_list 

        ##############################
        #  Make constraints
        ##############################

        print_verbose('make constraints',verbose)

        # TODO: renumber  these with the final numbering

        # note that the observations show up within model.act_windids as well, so we also constraint route scheduled DV by the real available DV from each observation
        def c1_rule( model,a):
            return (model.par_act_capacity[a]*model.var_activity_utilization[a] -
                        sum(model.par_dmr_act_dv[p,a]*model.var_dmr_utilization[p] 
                            for p in model.par_dmr_subscrs_by_act[a]) 
                    >= 0)
        model.c1 =pe.Constraint ( model.act_windids,  rule=c1_rule)

        # this constraint forces all activity indicators along a route to be high if that route is picked. From the other perspective, if an activity is not picked, then all route indicators through it must be low (not required, but the MILP branch and cut algorithm can take advantage of this to search through binary variables more efficiently)
        # it's not entirely clear to me though how helpful this constraint is - on a 30 sat model with 100 obs targets and 7 GS (1879 routes input to act sched) things seems to run slower by fractions of a second with this constaint (takes about 10 seconds to solve). Probably more helpful for larger models though...
        model.c1b  = pe.ConstraintList()
        for a in model.act_windids:
            for p in model.par_dmr_subscrs_by_act[a]:
                model.c1b.add(model.var_act_indic[a] >= model.var_dmr_indic[p]) 

        def c2_rule( model,p):
            return model.par_dmr_dv[p]*model.var_dmr_utilization[p] >= model.par_min_obs_dv_dlnk_req*model.var_dmr_indic[p]
        model.c2 =pe.Constraint ( model.dmr_ids,  rule=c2_rule)

        def c3_rule( model,a):
            return model.var_act_indic[a] >=  model.var_activity_utilization[a]
        model.c3 =pe.Constraint ( model.act_windids,  rule=c3_rule)  

        def c3c_rule( model,p):
            return model.var_dmr_indic[p] >=  model.var_dmr_utilization[p]
        model.c3c =pe.Constraint ( model.dmr_ids,  rule=c3c_rule)

        print_verbose('make overlap constraints',verbose)

        #  intra-satellite activity overlap constraints [4],[5],[5b]
        #  well, 5B is activity minimum time duration
        model.c4_5  = pe.ConstraintList() # this now contains all of the activity overlap constraints
        model.c5b  = pe.ConstraintList()
        # pass the model objects getter function so it can be called in place
        (self.c5b_binding_exprs_by_act,
            self.c4_5_binding_exprs_by_act) =  self.gen_intra_sat_act_overlap_constraints(
                model.c4_5,
                model.c5b,
                sats_acts,
                self.get_act_model_objs,
                constraint_violation_model_objs
            )

        # inter-satellite downlink overlap constraints [9],[10]
        model.c9_10  = pe.ConstraintList()
        self.c9_10_binding_exprs_by_act = self.gen_inter_sat_act_overlap_constraints(
            model.c9_10,sats_dlnks,
            self.get_act_model_objs,
            constraint_violation_model_objs
        )

        print_verbose('make energy, data constraints',verbose)

        #  energy constraints [6]
        # todo: maybe this ought to be moved to the super class, but i don't anticipate this code changing much any time soon, so i'll punt that.
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
                                act_windid = all_act_windids_by_obj[act]
                                act_code = act.get_code(sat_indx)
                                activity_delta_e += (
                                    model.par_sats_edot_by_mode[sat_indx][act_code] 
                                    * model.var_activity_utilization[act_windid]
                                    * model.par_resource_delta_t
                                )

                            #  if the satellite is not in sunlight then we can't charge
                            elif type(act) == EclipseWindow:
                                charging = False

                    # add in charging energy contribution ( if possible)
                    # assume charging is constant in sunlight
                    charging_delta_e = model.par_sats_edot_by_mode[sat_indx]['orbit_insunlight_average_charging']*model.par_resource_delta_t if charging else 0

                    #  base-level satellite energy usage (not including additional activities)
                    base_delta_e = model.par_sats_edot_by_mode[sat_indx]['base']*model.par_resource_delta_t

                    # maximum bound of energy at current time step based on last time step
                    model.c6.add( model.var_sats_estore[sat_indx,tp_indx] <= 
                        model.var_sats_estore[sat_indx,tp_indx-1]
                        + activity_delta_e
                        + charging_delta_e
                        + base_delta_e
                    )

                    # minimum bound of energy at current time step based on last time step
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

                    # todo: this may be too slow to use == below. Change it back to >= and use a better approach to extract real data storage values in extract_resource_usage below? Leaving == for now because it means I can use these vars to extract output data usage values
                    # constrain the minimum data storage at this tp by the amount of data volume being buffered by each route passing through the sat
                    if routes_storing:
                        model.c7.add( model.var_sats_dstore[sat_indx,tp_indx] == sum(model.par_dmr_dv[p]*model.var_dmr_utilization[p] for p in routes_storing))
                    else:
                        model.c7.add( model.var_sats_dstore[sat_indx,tp_indx] == 0)



        #  observation latency score factor constraints [8]
        model.c8  = pe.ConstraintList()
        for o in model.obs_windids:
            dmrs_obs = model.par_dmr_subscrs_by_obs_act[o]

            #  sort the latency score factors for all the dmrs for this observation in increasing order -  important for constraint construction
            dmrs_obs.sort(key= lambda p: latency_sf_by_dmr_id[p])

            num_dmrs_obs = len(dmrs_obs)
            #  initial constraint -  score factor for this observation will be equal to zero if no dmrs for this obs were chosen
            model.c8.add( model.var_latency_sf_obs[o] <= 0 + self.big_M_lat * sum(model.var_dmr_indic[p] for p in dmrs_obs) )
            
            for dmr_obs_indx in range(num_dmrs_obs):
                #  add constraint that score factor for observation is less than or equal to the score factor for this dmr_obs_indx, plus any big M terms for any dmrs with larger score factors.
                #  what this does is effectively disable the constraint for the score factor for this dmr_obs_indx if any higher score factor dmrs were chosen 
                model.c8.add( model.var_latency_sf_obs[o] <= 
                    latency_sf_by_dmr_id[dmrs_obs[dmr_obs_indx]] + 
                    self.big_M_lat * sum(model.var_dmr_indic[p] for p in dmrs_obs[dmr_obs_indx+1:num_dmrs_obs]) )

                #  note: use model.c8[indx].expr.to_string()  to print out the constraint in a human readable form
                #                ^ USES BASE 1 INDEXING!!! WTF??
                

        # constrain utilization of existing routes that are within planning_fixed_end
        model.c11  = pe.ConstraintList()
        for p in model.dmr_ids:
            if p in fixed_routes_ids:
                # less than constraint because equality should be achievable (if we're only using existing routes that have all previously been scheduled and deconflicted together - which is the case for current version of GP), but want to allow route to lessen its utilization if a more valuable route is available. 
                model.c11.add( model.var_dmr_utilization[p] <= utilization_by_existing_route_id[p]) 
            
        print_verbose('make obj',verbose)


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
        rsrc_norm_f = len(decimated_tp_indcs) * len(model.sat_indcs)

        def obj_rule(model):
            total_dv_term = self.obj_weights['obs_dv'] * 1/model.par_total_obs_dv * sum(model.par_dmr_dv[p]*model.var_dmr_utilization[p] for p in model.dmr_ids) 

            latency_term = self.obj_weights['route_latency'] * 1/len(model.obs_windids) * sum(model.var_latency_sf_obs[o] for o in model.obs_windids)
            
            energy_margin_term = self.obj_weights['energy_storage'] * 1/rsrc_norm_f * sum(model.var_sats_estore[sat_indx,tp_indx]/model.par_sats_estore_max[sat_indx] for tp_indx in decimated_tp_indcs for sat_indx in model.sat_indcs)

            if len(min_var_inter_sat_act_constr_violation_list) > 0:
                inter_sat_act_constr_violations_term = self.obj_weights['inter_sat_act_constr_violations'] * 1/sum(min_var_inter_sat_act_constr_violation_list) * sum(model.var_inter_sat_act_constr_violations[indx] for indx in range(1,len(model.var_inter_sat_act_constr_violations)+1))
            else:
                inter_sat_act_constr_violations_term = 0

            if len(min_var_intra_sat_act_constr_violation_list) > 0:
                intra_sat_act_constr_violations_term = self.obj_weights['intra_sat_act_constr_violations'] * 1/sum(min_var_intra_sat_act_constr_violation_list) * sum(model.var_intra_sat_act_constr_violations[indx] for indx in range(1,len(model.var_inter_sat_act_constr_violations)+1))
            else:
                intra_sat_act_constr_violations_term = 0

            # from circinus_tools import debug_tools
            # debug_tools.debug_breakpt()

            return total_dv_term + latency_term + energy_margin_term - inter_sat_act_constr_violations_term - intra_sat_act_constr_violations_term
            
        model.obj = pe.Objective( rule=obj_rule, sense=pe.maximize )

        self.model_constructed = True

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
                    print (" %d %0.3f  %s"%(index, val, self.all_acts_by_windid[index]))

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

            elif str (v) =='var_inter_sat_act_constr_violations': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" %d %0.3f   %s"%(index, val, self.inter_sat_act_constr_violation_acts_list[index-1]))

            elif str (v) =='var_intra_sat_act_constr_violations': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" %d %0.3f   %s"%(index, val, self.intra_sat_act_constr_violation_acts_list[index-1]))

    def extract_utilized_routes( self, verbose = False):
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
        for p in self.model.dmr_ids:
            if pe.value(self.model.var_dmr_indic[p]) >= 1.0 - self.binary_epsilon:
                # scheduled_route =  copy_choice (self.routes_filt[p]) 
                scheduled_route =  self.routes_by_dmr_id[p]
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
                act_indx = self.all_act_windids_by_obj[wind]
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

        t_vals = [[] for sat_indx in range ( self.num_sats)]
        e_vals = [[] for sat_indx in range ( self.num_sats)]

        # note that this extraction uses the energy variables from the optimization, which are currently not constrained to be exactly equal to the energy delta from t-1 to t; they are merely bounded by it. Todo: extract concrete values based on activity execution times
        # TODO: this code feels super inefficient somehow.  make it better?
        for sat_indx, sat in enumerate (self.model.sat_indcs):
            last_tp_indx = 0
            for tp_indx in self.model.es_timepoint_indcs:

                if (tp_indx - last_tp_indx) < decimation_factor:
                    continue
                else:
                    e_vals[sat_indx].append(pe.value(self.model.var_sats_estore[sat,tp_indx]))

                    t_vals[sat_indx].append(self.es_time_getter_dc.get_tp_from_tp_indx(tp_indx,out_units='minutes'))

        energy_usage['time_mins'] = t_vals
        energy_usage['e_sats'] = e_vals


        data_usage = {}

        t_vals = [[] for sat_indx in range ( self.num_sats)]
        d_vals = [[] for sat_indx in range ( self.num_sats)]

        # note that these d storage values are inexact in time, because they're determined conservatively from the original start and end of every window in each dmr. Window start/end times are adjusted in the actual scheduling, but we'll still say that a sat is holding the dmr's dv for the original extent of the window. That's why the data storage plot can look a bit off from the timing of the executed windows.
        # todo: if a more precise accounting of the dv on each satellite at each time is needed, then one should go through all windows and do bookkeeping on when exactly DV is moved from sat to sat.

        # TODO: this code feels super inefficient somehow.  make it better?
        for sat_indx, sat in enumerate (self.model.sat_indcs):
            last_tp_indx = 0
            for tp_indx in self.model.ds_timepoint_indcs:

                if (tp_indx - last_tp_indx) < decimation_factor:
                    continue
                else:
                    d_vals[sat_indx].append(pe.value(self.model.var_sats_dstore[sat,tp_indx]))

                    t_vals[sat_indx].append(self.ds_time_getter_dc.get_tp_from_tp_indx(tp_indx,out_units='minutes'))

        data_usage['time_mins'] = t_vals
        data_usage['d_sats'] = d_vals

        return  energy_usage, data_usage


    def extract_schedule_reasoning( self, all_routes, verbose = False):

        if verbose:
            print ('------------------------------')
            print ('extract_schedule_reasoning()')

        reasons = OrderedDict([
            ("a", "scheduled successfully"),
            ("b", "Observation window capacity fully utilized"),
            ("c", "cross-link window capacity fully utilized"),
            ("d", "downlink window capacity fully utilized"),
            ("e", "data route utilization below minimum (symptomatic)"),
            # ("f", "(non-conclusive) intra-satellite or inter-satellite activity overlap, or minimum activity time constraint violated"),
            ("g", "energy storage too low"),
            ("h", "data storage too high"),
            ("i", "intra-satellite activity overlap constraint would be violated"),
            ("j", "inter-satellite downlink overlap constraint would be violated"),
            ("k", "route filtered before scheduling"),
            ("z", "no reason found")
        ])

        reasons_by_route = {}

        # see https://en.wikipedia.org/wiki/Slack_variable - a constraint is binding if it's slack variable is zero
        #  this essentially tells us, in the case where an activity is involved in a constraint, that constraint is  binding, and the activity indicator for this activity is zero,  that the activity has been "forced out" of being allowed to execute because the execution of another activity constrained it
        #  the binding expression is essentially the constraint with any big M terms removed. 
        def expression_is_binding(be):
            return pe.value(be) < 0.001  # i.e., basically zero

        for dmr in all_routes:

            dmr_id = dmr.ID
            reasons_by_route[dmr] = set()

            if dmr_id not in self.routes_filt:
                reasons_by_route[dmr].add('k')
                continue

            #  check if successfully scheduled
            if pe.value(self.model.var_dmr_indic[dmr_id]) >= 1.0 - self.binary_epsilon:
                reasons_by_route[dmr].add('a')

            else:

                #  check if data route data volume requirement is too low
                if pe.value(self.model.var_dmr_utilization[dmr_id]) * dmr.data_vol < self.min_obs_dv_dlnk_req: 
                    reasons_by_route[dmr].add('e')

                for act in dmr.get_winds():
                    #  check if window capacity is fully or very near fully utilized
                    if act.scheduled_data_vol >= act.data_vol - self.dv_epsilon:
                        if type(act) == ObsWindow:
                            reasons_by_route[dmr].add('b')
                        if type(act) == DlnkWindow:
                            reasons_by_route[dmr].add('c')
                        if type(act) == XlnkWindow:
                            reasons_by_route[dmr].add('d')

                    #  if activity indicator is not set high, this means the activity was not used at all.  in this case it's possible that a constraint which relies on the activity indicator was violated
                    if pe.value(self.model.var_act_indic[act.window_ID]) <= 1.0 - self.binary_epsilon:
                        # reasons_by_route[dmr].add('f')
                        if act in self.c4_5_binding_exprs_by_act.keys():
                            binds = self.c4_5_binding_exprs_by_act[act]
                            for bind_expr in binds:
                                if expression_is_binding(bind_expr):
                                    #  constraint is binding and this activity is not being executed, so it should be the case
                                    reasons_by_route[dmr].add('i')

                        if act in self.c9_10_binding_exprs_by_act.keys():
                            binds = self.c9_10_binding_exprs_by_act[act]
                            for bind_expr in binds:
                                if expression_is_binding(bind_expr):
                                    #  constraint is binding and this activity is not being executed, so it should be the case
                                    reasons_by_route[dmr].add('j')


                        # note: can't really gain any information from c5b_binding_exprs_by_act

                    def get_act_resource_usage(act,model_resource,sat_indx):
                        tp_indx_center = self.ds_time_getter_dc.get_tp_indx_pre_t(act.center)
                        return pe.value(model_resource[sat_indx,tp_indx_center])


                    # check if energy storage is too low at an activity
                    estore_too_low_factor = 1.05 # assume if we're within certain % of energy lower bound then the activity is too constrained
                    if type(act) == ObsWindow or type(act) == DlnkWindow:
                        sat_estore = get_act_resource_usage(act,self.model.var_sats_estore,act.sat_indx)
                        if sat_estore < estore_too_low_factor * self.model.par_sats_estore_min[act.sat_indx]:
                            reasons_by_route[dmr].add('g')
                    if type(act) == XlnkWindow:
                        sat_estore = get_act_resource_usage(act,self.model.var_sats_estore,act.sat_indx)
                        if sat_estore < estore_too_low_factor * self.model.par_sats_estore_min[act.sat_indx]:
                            reasons_by_route[dmr].add('g')
                        xsat_estore = get_act_resource_usage(act,self.model.var_sats_estore,act.xsat_indx)
                        if xsat_estore < estore_too_low_factor * self.model.par_sats_estore_min[act.xsat_indx]:
                            reasons_by_route[dmr].add('g')


                    # check if data storage is too high at an activity
                    # if data storage is within the minimum data route data volume requirement of the maximum, then it's too high
                    if type(act) == ObsWindow or type(act) == DlnkWindow:
                        sat_dstore = get_act_resource_usage(act,self.model.var_sats_dstore,act.sat_indx)
                        if sat_dstore > self.model.par_sats_dstore_max[act.sat_indx] - self.min_obs_dv_dlnk_req:
                            reasons_by_route[dmr].add('h')
                    if type(act) == XlnkWindow:
                        sat_dstore = get_act_resource_usage(act,self.model.var_sats_dstore,act.sat_indx)
                        if sat_dstore > self.model.par_sats_dstore_max[act.sat_indx] - self.min_obs_dv_dlnk_req:
                            reasons_by_route[dmr].add('h')
                        xsat_dstore = get_act_resource_usage(act,self.model.var_sats_dstore,act.xsat_indx)
                        if xsat_dstore > self.model.par_sats_dstore_max[act.xsat_indx] - self.min_obs_dv_dlnk_req:
                            reasons_by_route[dmr].add('h')

 
            #  if we didn't find a reason, mark it. ( this is not great)
            if len(reasons_by_route[dmr]) == 0:
                reasons_by_route[dmr].add('z')


        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()


        all_reasons = []
        for dr_reasons in reasons_by_route.values():
            all_reasons+= list(dr_reasons)

        count_by_reason = OrderedDict()
        for reason in reasons.keys():
            count_by_reason[reason] = all_reasons.count(reason)

        num_successful = count_by_reason['a']
        num_unsuccessful = len(all_routes) - num_successful
        if verbose:
            print('num_routes: %d'%(len(all_routes)))
            print('num_successful: %d'%(num_successful))
            print('num_unsuccessful: %d'%(num_unsuccessful))
            print('#\t  reason')
            for reason_code,count in count_by_reason.items():
                print('%d\t: %s'%(count,reasons[reason_code]))


        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        return reasons_by_route











