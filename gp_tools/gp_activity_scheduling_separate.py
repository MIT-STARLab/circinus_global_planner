# contains model and Solver for global planner activity scheduling capability
# 
# @author Kit Kennedy
#

from  datetime import timedelta
from copy import  deepcopy
from math import ceil
from collections import OrderedDict,Counter

from pyomo import environ  as pe

import numpy as np

from circinus_tools  import time_tools as tt
from circinus_tools  import  constants as const
from circinus_tools.scheduling.custom_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow
from circinus_tools.scheduling.schedule_objects import Dancecard
from circinus_tools.scheduling.routing_objects import DataMultiRoute
from .gp_activity_scheduling_super import  GPActivityScheduling
from . import gp_general_tools as gp_gen

from circinus_tools import debug_tools

def print_verbose(string,verbose=False):
    if verbose:
        print(string)

class GPActivitySchedulingSeparate(GPActivityScheduling):
    """Non-optimal GP activity scheduling, to be run after performing route selection"""
    
    def __init__(self,gp_params):
        """initializes based on parameters
        
        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """

        super().__init__(gp_params)

        # number of routes found after filtering
        self.num_routes_filt = const.UNASSIGNED

        self.route_utilization_epsilon = 0.001  # 0.1% of route's DV

    def get_stats(self,verbose=True):

        num_winds_per_route = [len(dmr.get_winds()) for dmr in self.routes_filt]
        num_routes_by_act = {self.all_acts_by_windid[act_windid]:len(self.dmr_ids_by_act_windid[act_windid]) for act_windid in self.all_acts_windids}

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

    @staticmethod
    def get_act_model_objs(act,model):
        """ get the pyomo model objects used for modeling activity utilization."""
        # note: can be overridden in subclass - this function provides an example implementation

        model_objs_act = {
            'act_object': act,
            'var_dv_utilization': model.var_activity_utilization[act.window_ID]*model.par_act_capacity[act.window_ID],
            'par_dv_capacity': model.par_act_capacity[act.window_ID],
            'var_act_indic': model.var_act_indic[act.window_ID],
        }

        return model_objs_act

    def get_activity_structs( self,routes_filt):

        #  all activities are uniquely indexed. these structures keep track of those, and the mapping to activity objects
        all_acts_windids = set()
        mutable_acts_windids = set()
        # acts that are in planning window (both mutable and fixed acts)
        planwind_acts_windids = set()
        fixed_acts_windids = set()

        sats_planwind_acts = [[] for sat_indx in range (self.num_sats)]
        sats_planwind_dlnks = [[] for sat_indx in range (self.num_sats)]

        #  note that the rest of the structures consider both fixed and mutable activities

        #  these structures are for lookup in both directions
        all_acts_by_windid = {}

        # these structures keep track of the subset of unique indices that correspond to observations and links. Also we keep track of what data routes correspond to an activity
        all_lnk_windids = set()
        all_obs_windids = set()
        capacity_by_link_act_windid = {}
        capacity_by_obs_act_windid = {}
        capacity_by_act_windid = {}

        # now dmr structs
        # dmr is data multi-route
        dmr_ids_by_link_act_windid = {}
        dmr_ids_by_obs_act_windid = {}
        dmr_ids_by_act_windid = {}

        for dmr in routes_filt:
            for act in dmr.get_winds():

                # We have already filtered routes, but we also need to filter activities because there may be activities from existing routes that are outside of the planning window. note that , we do want to include acts outside of plan wind for other calculations in the model
                act_in_plan_wind = gp_gen.wind_in_planning_window(self,act,plan_wind_opt='whole')
                # also filter for acts that can be modified. We do not want to enforce constraints on those windows for which utilization cannot change
                act_is_mutable = gp_gen.wind_in_planning_window(self,act,plan_wind_opt='mutable')

                act_windid = act.window_ID

                # if we haven't yet seen this activity, then add it to bookkeeping
                #  note that this also filters out any activity object copies ( which is important for imposing activity scheduling constraints)
                if not act_windid in all_acts_windids:
                    # use activity window ID as a unique id for window in the model
                    all_acts_windids.add(act_windid)
                    if act_is_mutable: mutable_acts_windids.add(act_windid)
                    if act_in_plan_wind: planwind_acts_windids.add(act_windid)
                    if (act_in_plan_wind and not act_is_mutable): fixed_acts_windids.add(act_windid)
                    all_acts_by_windid[act_windid] = act
                    dmr_ids_by_act_windid[act_windid] = []
                    dmr_ids_by_act_windid[act_windid].append (dmr.ID)
                    capacity_by_act_windid[act_windid] = act.original_data_vol

                    # also need to add it to the list and dictionary for observations
                    if type(act) == ObsWindow:
                        if act_in_plan_wind: sats_planwind_acts[act.sat_indx].append(act)
                        all_obs_windids.add(act_windid)
                        dmr_ids_by_obs_act_windid[act_windid] = []
                        dmr_ids_by_obs_act_windid[act_windid].append (dmr.ID)
                        capacity_by_obs_act_windid[act_windid] = act.original_data_vol

                    # also need to add it to the list and dictionary for links
                    if type(act) == DlnkWindow:
                        all_lnk_windids.add(act_windid)
                        dmr_ids_by_link_act_windid[act_windid] = []
                        dmr_ids_by_link_act_windid[act_windid].append (dmr.ID)
                        capacity_by_link_act_windid[act_windid] = act.original_data_vol
                        # grab the dlnks for each sat too, while we're looping through
                        if act_in_plan_wind: sats_planwind_dlnks[act.sat_indx].append(act)
                        if act_in_plan_wind: sats_planwind_acts[act.sat_indx].append(act)

                    if type(act) == XlnkWindow:
                        all_lnk_windids.add(act_windid)
                        dmr_ids_by_link_act_windid[act_windid] = []
                        dmr_ids_by_link_act_windid[act_windid].append (dmr.ID)
                        capacity_by_link_act_windid[act_windid] = act.original_data_vol
                        if act_in_plan_wind: sats_planwind_acts[act.sat_indx].append(act)
                        if act_in_plan_wind: sats_planwind_acts[act.xsat_indx].append(act)

                #  if we have already seen the activity,  then just need to update the appropriate structures
                else:
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
            sats_planwind_acts[sat_indx].sort(key=lambda x: x.center)
            sats_planwind_dlnks[sat_indx].sort(key=lambda x: x.center)

        return sats_planwind_acts,sats_planwind_dlnks,all_acts_windids,planwind_acts_windids,mutable_acts_windids,fixed_acts_windids,dmr_ids_by_act_windid,capacity_by_act_windid,all_acts_by_windid,all_obs_windids,dmr_ids_by_obs_act_windid,capacity_by_obs_act_windid,all_lnk_windids,dmr_ids_by_link_act_windid,capacity_by_link_act_windid
                    

    def filter_routes( self,selected_routes,existing_routes):

        routes_filt = []
        existing_routes_before_planning_fixed_end = []

        def route_has_too_short_act(rt):
            for wind in rt.get_winds():
                if wind.duration.total_seconds() < self.act_timing_helper.get_act_min_duration(wind):
                    return True

        def process_dmr(dmr,plan_wind_opt,filter_opt):
            dmr_start = dmr.get_start(time_opt='original')
            dmr_end = dmr.get_end(time_opt='original')

            # check if all act windows in route are completely within the planning window. Pass if not.
            if filter_opt=='totally_within' and not gp_gen.dr_in_planning_window(self,dmr,plan_wind_opt): 
                pass
            # check if at least one act window in route is partially within the planning window. Pass if not.
            elif filter_opt=='partially_within' and (dmr_end < self.planning_start_dt or dmr_start > self.planning_end_dt):
                pass
            elif route_has_too_short_act(dmr):
                # print('discarding too short dlnk window')
                pass
            else:
                routes_filt.append (dmr)
                return True

            return False
            
        # want new routes to be entirely within non-fixed planning window (otherwise we could miss enforcing all of the constraints imposed by the route)
        for dmr in selected_routes:
            process_dmr(dmr,'mutable',"totally_within")

        # want existing routes to be at least partially within whole planning window, so that we will be able to later filter for all act windows on which they impose constraints (but note that we don't have to worry anymore about the constraints on the route, because the input utilization number for the route accounts for it) 
        for dmr in existing_routes:
            added = process_dmr(dmr,'whole',"partially_within")

            # mark those  existing routes that start before the fixed planning window end.
            dmr_start = dmr.get_start(time_opt='original')
            if added and dmr_start <= self.planning_fixed_end_dt:
                existing_routes_before_planning_fixed_end.append(dmr)

        return routes_filt, existing_routes_before_planning_fixed_end

    def make_model ( self,selected_routes, existing_routes, utilization_by_existing_route_id, ecl_winds, verbose = True):
        # important assumption: all activity window IDs are unique!

        #  first run pre-check on satellite states, to make sure that we can actually solve the MILP
        self.run_sat_state_precheck(ecl_winds)
        print_verbose('Passed state pre-check',verbose)

        model = pe.ConcreteModel()
        self.model = model

        ##############################
        # filtering/checking of selected routes
        ##############################

        # verify that each route is unique. If duplicate routes are present, would cause problems in constraints (e.g. a route would require double its actual needed throuput in a given act)
        existing_routes_set = set(existing_routes)
        #  make sure no existing routes are duplicated
        assert(len(existing_routes) == len(existing_routes_set))
        duplicated_new_routes = [dmr for dmr, count in Counter(selected_routes).items() if count > 1 and dmr not in existing_routes_set]
        #  this makes sure that for every route that is not in existing Routes (i.e., a new route), there is no duplicate for it
        assert( len(duplicated_new_routes) == 0)

        #  in the  route selection stage, we could have added some of the existing routes as "newly selected routes".  this is okay.  however, we do want to remove them from new routes right now in order to avoid constraint problems as mentioned above.
        new_routes = []
        #  count the number of times we removed an existing route from the new routes.  ideally all existing routes would be included in new routes (note: new routes are made at route selection step two).
        existing_routes_selected = set()
        for dmr in selected_routes:
            if dmr in existing_routes_set:
                existing_routes_selected.add(dmr)
            else:
                new_routes.append(dmr)

        # filter the routes to make sure that  none of their activities fall outside the scheduling window
        routes_filt,existing_routes_fixed = self.filter_routes(new_routes,existing_routes)
        routes_by_dmr_id = {dmr.ID:dmr for dmr in routes_filt}

        print_verbose('considering %d routes'%(len(routes_filt)),verbose)
        print_verbose('Fraction of existing routes included at RS step two: %d/%d'%(len(existing_routes_selected),len(existing_routes)),verbose)

        self.routes_filt = routes_filt
        self.existing_routes = existing_routes
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

        fixed_routes_ids = [dmr.ID for dmr in existing_routes_fixed]
        existing_routes_ids = [dmr.ID for dmr in existing_routes]
        sum_existing_route_utilization = sum(utilization_by_existing_route_id[rt_id] for rt_id in existing_routes_ids)
        self.fixed_routes_ids = fixed_routes_ids
        self.existing_routes_ids = existing_routes_ids

        ##############################
        #  Make indices/ subscripts
        ##############################

        try:
            (sats_planwind_acts,
                sats_planwind_dlnks,                
                all_acts_windids,
                planwind_acts_windids,
                mutable_acts_windids,
                fixed_acts_windids,
                dmr_ids_by_act_windid,
                capacity_by_act_windid,
                all_acts_by_windid,
                all_obs_windids,
                dmr_ids_by_obs_act_windid,
                capacity_by_obs_act_windid,
                all_lnk_windids,
                dmr_ids_by_link_act_windid,
                capacity_by_link_act_windid) =  self.get_activity_structs(routes_filt)

            latency_sf_by_dmr_id =  self.get_route_latency_score_factors(
                routes_by_dmr_id,
                dmr_ids_by_obs_act_windid
            )


            # debug_tools.debug_breakpt()

            fixed_wind_utilization_by_wind_id = {}
            for wind_id in fixed_acts_windids:    
                wind = all_acts_by_windid[wind_id]

                fixed_wind_utilization_by_wind_id[wind_id] = sum(rt.data_vol_for_wind(wind)*utilization_by_existing_route_id[rt.ID] for rt in existing_routes_fixed if wind in rt) / wind.original_data_vol

            # # verify that all acts found are within the planning window, otherwise we may end up with strange results
            # for sat_acts in sats_planwind_acts:
            #     for act in sat_acts:
            #         # note that these should not be original start/end. It's possible for downlinks that have an original start before the planning window to get added, in route selection, if they can still deliver data volume for an obs. 
            #         if act.start < self.planning_start_dt or act.end > self.planning_end_dt:
            #             raise RuntimeWarning('Activity is out of planning window range (start %s, end %s): %s'%(self.planning_start_dt,self.planning_end_dt,act))

            # construct a set of dance cards for every satellite, 
            # each of which keeps track of all of the activities of satellite 
            # can possibly execute at any given time slice delta T. 
            # this is for constructing energy storage constraints
            # using resource_delta_t_s because this dancecard is solely for use in constructing resource constraints
            # note that these dancecards will baloon in size pretty quickly as the planning_end_dt increases. However most of the complexity they introduce is before planning_end_obs,xlnk_dt because that's the horizon where obs,xlnk actitivities are included. After that there should only be sparse downlinks
            es_act_dancecards = [Dancecard(self.planning_start_dt,self.planning_end_dt,self.resource_delta_t_s,item_init=None,mode='timestep') for sat_indx in range (self.num_sats)]
            
            #  add windows to dance card, silenty dropping any time steps that appear outside of the planning window bounds ( we don't need to enforce resource constraints on out-of-bounds activities)
            def wind_time_getter_orig(wind,time_opt):
                if time_opt == 'start': return wind.original_start
                if time_opt == 'end': return wind.original_end
            def wind_time_getter_reg(wind,time_opt):
                if time_opt == 'start': return wind.start
                if time_opt == 'end': return wind.end

            for sat_indx in range (self.num_sats): 
                # note that we need to consider all acts in the plan window, not just mutables
                # note that we want to use original start/end times here because start/end might have been changed on a previous global planner run, if this is running in the const simulation context.
                es_act_dancecards[sat_indx].add_winds_to_dancecard(sats_planwind_acts[sat_indx],wind_time_getter_orig,drop_out_of_bounds=True)
                es_act_dancecards[sat_indx].add_winds_to_dancecard(ecl_winds[sat_indx],wind_time_getter_reg,drop_out_of_bounds=True)

            # this is for data storage
            # for each sat/timepoint, we store a list of all those data multi routes that are storing data on the sat at that timepoint
            ds_route_dancecards = [Dancecard(self.planning_start_dt,self.planning_end_dt,self.resource_delta_t_s,item_init=None,mode='timepoint') for sat_indx in range (self.num_sats)]
            
            # add data routes to the dancecard
            for dmr in routes_filt:
                # list of type routing_objects.SatStorageInterval
                dmr_ds_intervs = dmr.get_data_storage_intervals()

                for interv in dmr_ds_intervs:
                    # store the dmr object at this timepoint, silenty dropping any time points that appear outside of the planning window bounds ( we don't need to enforce resource constraints on out-of-bounds intervals)
                    ds_route_dancecards[interv.sat_indx].add_item_in_interval(dmr.ID,interv.start,interv.end,drop_out_of_bounds=True)

        except IndexError:
            raise RuntimeWarning('sat_indx out of range. Are you sure all of your input files are consistent? (including pickles)')        
        self.dmr_ids_by_act_windid = dmr_ids_by_act_windid
        self.all_acts_windids = all_acts_windids
        self.obs_windids = all_obs_windids
        self.lnk_windids = all_lnk_windids
        self.all_acts_by_windid = all_acts_by_windid
        self.mutable_acts_windids = mutable_acts_windids
        self.planwind_acts_windids = planwind_acts_windids

        #  subscript for each dmr (data multi route) p  (p index is a hold-over from when I called them paths)
        model.dmr_ids = pe.Set(initialize= routes_by_dmr_id.keys())
        #  these are all the routes that have acts occurring before the fixed planning window end.  we assume that these routes are no longer as "malleable" as other routes -  they have an upper limit on their utilization based on their existing utilization. this is because if any of these activities has already been executed, the route is now constrained by the throughput used from that act
        model.fixed_dmr_ids = pe.Set(initialize= fixed_routes_ids)
        model.existing_dmr_ids = pe.Set(initialize= existing_routes_ids)
        #  subscript for each activity a
        model.all_act_windids = pe.Set(initialize= all_acts_windids)
        #  subscript for each mutable activity a ( we can still change the activity's utilization)
        model.planwind_acts_windids = pe.Set(initialize= planwind_acts_windids)
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
        elif self.solver_name == 'glpk':
            # raise an error, because it could be misleading if someone changes the int feas tol in the inputs...
            raise NotImplementedError('glpk runs, but I have not yet figured out setting integer_feasibility_tolerance')
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
        model.par_obs_capacity = pe.Param(model.obs_windids,initialize =capacity_by_obs_act_windid)
        model.par_total_obs_capacity = sum(capacity_by_obs_act_windid.values())
        model.par_link_capacity = pe.Param(model.lnk_windids,initialize =capacity_by_link_act_windid)
        model.par_act_capacity = pe.Param(model.all_act_windids,initialize =capacity_by_act_windid)
        #  data volume for each data multi-route
        model.par_dmr_dv = pe.Param(model.dmr_ids,initialize ={ dmr.ID: dmr.data_vol for dmr in routes_filt})
        model.par_fixed_wind_utilization_by_wind_id = pe.Param(model.fixed_acts_windids,initialize ={ wind_id: util for wind_id,util in fixed_wind_utilization_by_wind_id.items()})
        #  data volume for each activity in each data multi-route

        #  stored data volume for each activity used by each data route. only store this for the activity if the activity was found to be within the planning window filter (hence, the act in self.all_act_windids_by_obj.keys() check)
        model.par_dmr_act_dv = pe.Param(
            model.dmr_ids,
            model.all_act_windids,
            initialize = { (dmr.ID,act.window_ID): 
                dmr.data_vol_for_wind(act) for dmr in routes_filt for act in dmr.get_winds() if act.window_ID in all_acts_windids 
            }
        )

        # each of these is essentially a dictionary indexed by link or obs act indx, with  the value being a list of dmr indices that are included within that act
        # these are valid indices into model.dmr_ids
        # model.par_dmr_subscrs_by_link_act = pe.Param(model.lnk_windids,initialize =dmr_ids_by_link_act_windid)
        model.par_dmr_subscrs_by_obs_act = pe.Param(model.obs_windids,initialize =dmr_ids_by_obs_act_windid)
        model.par_dmr_subscrs_by_act = pe.Param(model.all_act_windids,initialize =dmr_ids_by_act_windid)


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
        model.var_activity_utilization  = pe.Var (model.planwind_acts_windids, bounds =(0,1))
        # dmr utilization variable indicating how much of a dmr's capacity is used [2]
        model.var_dmr_utilization  = pe.Var (model.dmr_ids, bounds =(0,1))
        #  indicator variables for whether or not dmrs [3] and activities [4] have been chosen
        model.var_dmr_indic  = pe.Var (model.dmr_ids, within = pe.Binary)
        model.var_act_indic  = pe.Var (model.planwind_acts_windids, within = pe.Binary)

        # a utilization number for existing routes that will be bounded by the input existing route utilization (can't get more  "existing route" reward for a route than the route's previous utilization) [8]
        model.var_existing_dmr_utilization_reward  = pe.Var (model.existing_dmr_ids, bounds =(0,1))
        
        # satellite energy storage [5]
        model.var_sats_estore  = pe.Var (model.sat_indcs,  model.es_timepoint_indcs,  within = pe.NonNegativeReals)

        # satellite data storage (data buffers) [6]
        model.var_sats_dstore  = pe.Var (model.sat_indcs,  model.ds_timepoint_indcs,  within = pe.NonNegativeReals)

        # Reward factor for latency for a given observation [7]
        model.var_latency_sf_obs = pe.Var (model.obs_windids,  bounds = (0,1.0))
        
        allow_act_timing_constr_violations = False
        if allow_act_timing_constr_violations:
            print('allow_act_timing_constr_violations is True')

        #  variables for handling the allowance of inter-activity timing constraint violations. these are only generated if allow_act_timing_constr_violations is True
        model.var_intra_sat_act_constr_violations = pe.VarList()
        model.var_inter_sat_act_constr_violations = pe.VarList()
        model.intra_sat_act_constr_bounds  = pe.ConstraintList()
        model.inter_sat_act_constr_bounds  = pe.ConstraintList()

        self.constraint_violation_model_objs['var_intra_sat_act_constr_violations'] = model.var_intra_sat_act_constr_violations
        self.constraint_violation_model_objs['var_inter_sat_act_constr_violations'] = model.var_inter_sat_act_constr_violations
        self.constraint_violation_model_objs['intra_sat_act_constr_bounds'] = model.intra_sat_act_constr_bounds
        self.constraint_violation_model_objs['inter_sat_act_constr_bounds'] = model.inter_sat_act_constr_bounds

        ##############################
        #  Make constraints
        ##############################

        print_verbose('make constraints',verbose)

        # TODO: renumber  these with the final numbering

        # note that the observations show up within model.all_act_windids as well, so we also constraint route scheduled DV by the real available DV from each observation
        def c1_rule( model,a):
            return (model.par_act_capacity[a]*model.var_activity_utilization[a] -
                        sum(model.par_dmr_act_dv[p,a]*model.var_dmr_utilization[p] 
                            for p in model.par_dmr_subscrs_by_act[a]) 
                    >= 0)
        model.c1 =pe.Constraint ( model.planwind_acts_windids,  rule=c1_rule)

        # this constraint forces all activity indicators along a route to be high if that route is picked. From the other perspective, if an activity is not picked, then all route indicators through it must be low (not required, but the MILP branch and cut algorithm can take advantage of this to search through binary variables more efficiently)
        # it's not entirely clear to me though how helpful this constraint is - on a 30 sat model with 100 obs targets and 7 GS (1879 routes input to act sched) things seems to run slower by fractions of a second with this constaint (takes about 10 seconds to solve). Probably more helpful for larger models though...
        model.c1b  = pe.ConstraintList()
        for a in model.planwind_acts_windids:
            for p in model.par_dmr_subscrs_by_act[a]:
                model.c1b.add(model.var_act_indic[a] >= model.var_dmr_indic[p]) 

        def c2_rule( model,p):
            return model.par_dmr_dv[p]*model.var_dmr_utilization[p] >= model.par_min_obs_dv_dlnk_req*model.var_dmr_indic[p]
        model.c2 =pe.Constraint ( model.dmr_ids,  rule=c2_rule)

        def c3_rule( model,a):
            return model.var_act_indic[a] >=  model.var_activity_utilization[a]
        model.c3 =pe.Constraint ( model.planwind_acts_windids,  rule=c3_rule)  

        # todo: this probably should not be enforced...
        # def c3c_rule( model,p):
        #     return model.var_dmr_indic[p] >=  model.var_dmr_utilization[p]
        # model.c3c =pe.Constraint ( model.dmr_ids,  rule=c3c_rule)


        def c3d_rule( model,a):
            return model.par_fixed_wind_utilization_by_wind_id[a] ==  model.var_activity_utilization[a]
        model.c3d =pe.Constraint ( fixed_acts_windids,  rule=c3d_rule)  

        print_verbose('make overlap constraints',verbose)

        #  intra-satellite activity overlap constraints [4],[5]
        model.c4_5  = pe.ConstraintList()
        # pass the model objects getter function so it can be called in place
        self.c4_5_binding_exprs_by_act =  self.gen_intra_sat_act_overlap_constraints(
            model,
            model.c4_5,
            sats_planwind_acts,
            self.num_sats,
            self.get_act_model_objs
        )

        #  5B is activity minimum time duration
        model.c5b  = pe.ConstraintList()
        # pass the model objects getter function so it can be called in place
        self.c5b_binding_exprs_by_act =  self.gen_sat_act_duration_constraints(
            model,
            model.c5b,
            sats_planwind_acts,
            self.num_sats,
            self.get_act_model_objs
        )

        # inter-satellite downlink overlap constraints [9],[10]
        model.c9_10  = pe.ConstraintList()
        self.c9_10_binding_exprs_by_act = self.gen_inter_sat_act_overlap_constraints(
            model,
            model.c9_10,
            sats_planwind_dlnks,
            self.num_sats,
            self.get_act_model_objs
        )

        print_verbose('make energy, data constraints',verbose)


        #  energy constraints [6]
        # todo: maybe this ought to be moved to the super class, but i don't anticipate this code changing much any time soon, so i'll punt that.
        model.c6  = pe.ConstraintList()
        for sat_indx in range (self.num_sats): 

            charge_eff = self.sats_batt_charge_eff[sat_indx]
            discharge_eff = self.sats_batt_discharge_eff[sat_indx]
            discharge_factor = 1/discharge_eff

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
                                act_code = act.get_e_dot_codename(sat_indx)
                                activity_delta_e += (
                                    model.par_sats_edot_by_mode[sat_indx][act_code] 
                                    * model.var_activity_utilization[act.window_ID]
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
                    # assume 100% of energy can be delivered if we're not charging the battery
                    model.c6.add( model.var_sats_estore[sat_indx,tp_indx] <= 
                        model.var_sats_estore[sat_indx,tp_indx-1]
                        + discharge_factor*activity_delta_e
                        + charging_delta_e
                        + discharge_factor*base_delta_e
                    )

                    # assume max recharge amount limited by charging efficiency.
                    model.c6.add( model.var_sats_estore[sat_indx,tp_indx] <= 
                        model.var_sats_estore[sat_indx,tp_indx-1]
                        + charge_eff * charging_delta_e
                    )

                    # minimum bound of energy at current time step based on last time step
                    model.c6.add( model.var_sats_estore[sat_indx,tp_indx] >= 
                        model.var_sats_estore[sat_indx,tp_indx-1]
                        + discharge_factor * activity_delta_e
                        + discharge_factor * base_delta_e
                    )


        #  data storage constraints [7]
        # note that these constraints are overly conservative, because the data storage intervals calculated for each route are based on the original start and end times of activities, not on the start and end times adjusted for the activity's utilization. Unfortunately that's a difficulty of how we discretized the decision making here...
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
        #  note this is for all observation window IDs, not just mutable ones
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
        for p in model.fixed_dmr_ids:
            # less than constraint because equality should be achievable (if we're only using existing routes that have all previously been scheduled and deconflicted together - which is the case for current version of GP), but want to allow route to lessen its utilization if a more valuable route is available. 
            #  add in an epsilon at the end, because it may be that the utilization number was precisely chosen to meet the minimum data volume requirement -  don't want to not make the minimum data volume requirement this time because of round off error
            model.c11.add( model.var_dmr_utilization[p] <= utilization_by_existing_route_id[p] + self.fixed_utilization_epsilon) 
        
        # constrain utilization reward of existing routes by the input utilization numbers
        model.c12  = pe.ConstraintList()
        for p in model.existing_dmr_ids:
            # constrain the reward utilization (used in obj function) by the input utilization for this existing route
            model.c12.add( model.var_existing_dmr_utilization_reward[p] <= utilization_by_existing_route_id[p] ) 
            # also constrain by the current utilization in the optimization
            model.c12.add( model.var_existing_dmr_utilization_reward[p] <= model.var_dmr_utilization[p] ) 

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
            #  note the first two objectives are for all observations, not just mutable observations

            # obj [1]
            total_dv_term = self.obj_weights['obs_dv'] * 1/model.par_total_obs_capacity * sum(model.par_dmr_dv[p]*model.var_dmr_utilization[p] for p in model.dmr_ids) 

            # obj [2]
            latency_term = self.obj_weights['route_latency'] * 1/len(model.obs_windids) * sum(model.var_latency_sf_obs[o] for o in model.obs_windids)
            
            # obj [5]
            energy_margin_term = self.obj_weights['energy_storage'] * 1/rsrc_norm_f * sum(model.var_sats_estore[sat_indx,tp_indx]/model.par_sats_estore_max[sat_indx] for tp_indx in decimated_tp_indcs for sat_indx in model.sat_indcs)

            # obj [6]
            existing_routes_term = self.obj_weights['existing_routes'] * 1/sum_existing_route_utilization * sum(model.var_existing_dmr_utilization_reward[p] for p in model.existing_dmr_ids) if len(model.existing_dmr_ids) > 0 else 0

            if len(self.min_var_inter_sat_act_constr_violation_list) > 0:
                inter_sat_act_constr_violations_term = self.obj_weights['inter_sat_act_constr_violations'] * 1/sum(self.min_var_inter_sat_act_constr_violation_list) * sum(model.var_inter_sat_act_constr_violations[indx] for indx in range(1,len(model.var_inter_sat_act_constr_violations)+1))
            else:
                inter_sat_act_constr_violations_term = 0

            if len(self.min_var_intra_sat_act_constr_violation_list) > 0:
                intra_sat_act_constr_violations_term = self.obj_weights['intra_sat_act_constr_violations'] * 1/sum(self.min_var_intra_sat_act_constr_violation_list) * sum(model.var_intra_sat_act_constr_violations[indx] for indx in range(1,len(model.var_inter_sat_act_constr_violations)+1))
            else:
                intra_sat_act_constr_violations_term = 0

            # from circinus_tools import debug_tools
            # debug_tools.debug_breakpt()

            return total_dv_term + latency_term + energy_margin_term + existing_routes_term - inter_sat_act_constr_violations_term - intra_sat_act_constr_violations_term
            
        model.obj = pe.Objective( rule=obj_rule, sense=pe.maximize )

        self.model_constructed = True

        # print(self.planning_start_dt)
        # print(self.planning_end_dt)
        # print([rt.get_obs().end for rt in  selected_routes])
        # print('all minus mutable')
        # print([self.all_acts_by_windid[a] for a in (set(all_acts_windids)-set(mutable_acts_windids))])
        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

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
        #  note: this should be the only place where scheduled data volume attributes are updated

        # Keeps track of which routes were chosen to be executed
        scheduled_routes = []
        #  keeps track of all routes that were present in activity scheduling, both of those that ended up being chosen (scheduled), and those that ended up not being chosen
        all_updated_routes = []

        if verbose:
            print ('utilized routes:')


        # figure out which dmrs were used, and add the scheduled data volume to each
        for p in self.model.dmr_ids:
            route =  self.routes_by_dmr_id[p]
            route.set_scheduled_dv_frac(pe.value(self.model.var_dmr_utilization[p]))  #* self.model.par_dmr_dv[p])
            all_updated_routes. append (route)
            
            #  if it was actually scheduled; that is, the indicator is greater than zero
            if pe.value(self.model.var_dmr_utilization[p]) >= self.route_utilization_epsilon:
                scheduled_routes. append (route)

                if verbose:
                    print(route)


        scheduled_dv_by_wind = {}

        # examine the schedulable data volume for every activity window, checking as we go that the data volume is sufficient for at least the route in which the window is found
        #  note that this code is slightly inefficient because it might duplicate windows across routes. that's fine though, because we're thorough in checking across all routes
        # note: dmr is for DataMultiRoute
        wind_sched_dv_check = {}
        for dmr in all_updated_routes:
            # wind may get set multiple times due to Windows appearing across routes, but that's not really a big deal
            for wind in dmr.get_winds():
                scheduled_dv_by_wind[wind] = 0

                #  only update windows that are mutable
                if wind.window_ID in self.mutable_acts_windids:
                    # wind_sched_dv_check[wind] = wind.data_vol * pe.value(self.model.var_act_indic[act_indx])
                    wind_sched_dv_check[wind] = wind.original_data_vol * pe.value(self.model.var_activity_utilization[wind.window_ID])
                    #  initialize this while we're here
                    # wind.scheduled_data_vol = 0
                else:
                    #  if it's not a mutable window it should already have scheduled data volume assigned
                    # assert(wind.scheduled_data_vol != const.UNASSIGNED)
                    wind_sched_dv_check[wind] = 0

            #  update the underlying data routes for every data multi-route
            #  todo:  should this be included in set_scheduled_dv_frac() call above?
            # similar to winds, set dr data vols to 0
            for dr in dmr.data_routes:
                if not dr in self.existing_routes:
                    pass
                    # dmrs should not be sharing drs across themselves, so no scheduled dv should be seen yet
                    # todo: there seems to be a bug here...to debug.
                    # assert(dr.scheduled_dv == const.UNASSIGNED)
                dr.scheduled_dv = dmr.scheduled_dv_by_dr[dr] 


        #  now we want to mark the real scheduled data volume for every window. We need to do this separately because the model.var_act_indic continuous variables only give an upper bound on the data volume for an activity. we only actually need to use as much data volume as the data routes want to push through the window
        #  add data volume for every route passing through every window
        for dmr in all_updated_routes:
            for wind in dmr.get_winds():
                #  only update windows that are mutable
                if wind.window_ID in self.mutable_acts_windids:
                    scheduled_dv_by_wind[wind] += dmr.scheduled_dv_for_wind(wind)
                else:
                    wind_sched_dv_check[wind] +=  dmr.scheduled_dv_for_wind(wind)


        # update the window beginning and end times based upon their amount of scheduled data volume
        # keep track of which windows we've updated, because we should only update once
        updated_winds = set()
        for dmr in all_updated_routes:
            # validate the data multi route (and in turn, the scheduled data vols of all the data routes under it)
            dmr.validate(self.act_timing_helper)

            for wind in dmr.get_winds():

                if wind.window_ID in self.mutable_acts_windids:
                    #  this check should be at least as big as the scheduled data volume as calculated from all of the route data volumes. (it's not constrained from above, so it could be bigger)
                    if wind_sched_dv_check[wind] < scheduled_dv_by_wind[wind] - self.dv_epsilon:
                        raise RuntimeWarning('inconsistent activity scheduling results, data volumes mismatch, mutable window')


                # if the window is not mutable, then it only should be here because of an existing data route. in that case we know the previous utilization for the window.  check the current schedule datable volume for window against the previous utilization
                else:
                    previous_dv_utilization = sum(self.utilization_by_existing_route_id[dmr.ID]*dmr.data_vol_for_wind(wind) for dmr in self.existing_routes if wind in dmr.get_winds())

                    # check if somehow data routes through fixed window have tried to grab more capacity than was scheduled for it
                    if wind_sched_dv_check[wind] > previous_dv_utilization + self.dv_epsilon:
                        raise RuntimeWarning('inconsistent activity scheduling results, data volumes mismatch, fixed window. Previous dv scheduled: %f, current scheduled %f. Verify that fixed_utilization_epsilon (%f) is not too large relative to dv_epsilon (%f)'%(previous_dv_utilization,wind_sched_dv_check[wind],self.fixed_utilization_epsilon,self.dv_epsilon))

                #  also check that we're not scheduling too much data volume from the window ( check this after we already verified data volume usage relative to previous utilization, so we see that error first -  helps to separate out that specific case)
                if wind_sched_dv_check[wind] >= wind.original_data_vol + self.dv_epsilon:
                    raise RuntimeWarning('too much data volume was scheduled for window %s'%(wind))


                #  only update windows that are mutable
                if not wind.window_ID in self.mutable_acts_windids:
                    continue

                # if it hasn't been updated yet
                if not wind in updated_winds:
                    # note that the line below seems like it may break the scheduled times for activities by specifying a minimum activity duration. however, this minimum activity duration is already accounted for in scheduling constraints
                    wind.scheduled_data_vol = scheduled_dv_by_wind[wind]
                    wind.update_duration_from_scheduled_dv (min_duration_s=self.act_timing_helper.get_act_min_duration(wind))
                    updated_winds.add(wind)


        return scheduled_routes,all_updated_routes

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
                    #  only consider windows that were in the planning window
                    if not act in self.planwind_acts_windids:
                        continue

                    #  check if window capacity is fully or very near fully utilized
                    if act.scheduled_data_vol >= act.original_data_vol - self.dv_epsilon:
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
                        tp_indx_center_at_act = self.ds_time_getter_dc.get_tp_indx_pre_t(act.center)
                        return pe.value(model_resource[sat_indx,tp_indx_center])

                    # note that the energy resource usage check below is not very robust, because the effects of energy usage by an activity are not necessarily instantaneous - they can overconstrain the satellite at a later time. Case in point: satellite collects obs data for too long at time t = 10 mins, while in eclipse. It has enough energy to stay above min bound for the whole activity, and a while afterward. But at time t = 20 mins, right before the end of the eclipse, energy goes below lower bound. If we had executed the obs and spent a little less energy at t = 10 mins, we'd be able to stay in the green, but we didn't, so we're hosed. (note this is a hypothetical - the scheduler is only capable of solving the problem if it reduces the obs time and stays within energy constraints.)

                    # check if energy storage is too low at an activity
                    estore_too_low_factor = 1.05 # assume if we're within certain % of energy lower bound then the activity is too constrained
                    if type(act) == ObsWindow or type(act) == DlnkWindow:
                        sat_estore = get_resource_usage_at_act(act,self.model.var_sats_estore,act.sat_indx)
                        if sat_estore < estore_too_low_factor * self.model.par_sats_estore_min[act.sat_indx]:
                            reasons_by_route[dmr].add('g')
                    if type(act) == XlnkWindow:
                        sat_estore = get_resource_usage_at_act(act,self.model.var_sats_estore,act.sat_indx)
                        if sat_estore < estore_too_low_factor * self.model.par_sats_estore_min[act.sat_indx]:
                            reasons_by_route[dmr].add('g')
                        xsat_estore = get_resource_usage_at_act(act,self.model.var_sats_estore,act.xsat_indx)
                        if xsat_estore < estore_too_low_factor * self.model.par_sats_estore_min[act.xsat_indx]:
                            reasons_by_route[dmr].add('g')

                    # data storage check is robust though, because data storage effects are instananeous, unlinke energy (see above)

                    # check if data storage is too high at an activity
                    # if data storage is within the minimum data route data volume requirement of the maximum, then it's too high
                    if type(act) == ObsWindow or type(act) == DlnkWindow:
                        sat_dstore = get_resource_usage_at_act(act,self.model.var_sats_dstore,act.sat_indx)
                        if sat_dstore > self.model.par_sats_dstore_max[act.sat_indx] - self.min_obs_dv_dlnk_req:
                            reasons_by_route[dmr].add('h')
                    if type(act) == XlnkWindow:
                        sat_dstore = get_resource_usage_at_act(act,self.model.var_sats_dstore,act.sat_indx)
                        if sat_dstore > self.model.par_sats_dstore_max[act.sat_indx] - self.min_obs_dv_dlnk_req:
                            reasons_by_route[dmr].add('h')
                        xsat_dstore = get_resource_usage_at_act(act,self.model.var_sats_dstore,act.xsat_indx)
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











