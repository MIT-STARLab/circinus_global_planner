# contains model and Solver for global planner activity scheduling capability
# 
# @author Kit Kennedy
#

from  datetime import timedelta
from copy import  deepcopy
from math import ceil
from copy import copy

from pyomo import environ  as pe
from pyomo import opt  as po

import numpy as np

from .gp_activity_scheduling_super import  GPActivityScheduling
from .routing_objects import DataRoute, DataMultiRoute
from circinus_tools  import time_tools as tt
from .custom_activity_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow
from .routing_objects import DataMultiRoute
from .schedule_objects import Dancecard
from circinus_tools  import  constants as const

class GPActivitySchedulingCoupled(GPActivityScheduling):
    """docstring for GP activity scheduling"""
    
    def __init__(self,gp_params):
        """initializes based on parameters
        
        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """

        super().__init__(gp_params)

        # if two times are separated by this amount, we call them effectively the same time
        self.T_s_time_epsilon_td = timedelta(seconds=self.scenario_timestep_s/2)

    def get_stats(self,verbose=True):

        stats = {}
        stats['num_acts'] = len ( self.all_acts_windids)
        stats['num_obs'] = len ( self.all_obs_windids)
        stats['num_dlnks'] = len ( self.all_dlnk_windids)
        stats['num_xlnks'] = len ( self.all_xlnk_windids)
        stats['num_variables'] = self.model.nvariables ()
        stats['num_nobjectives'] = self.model.nobjectives ()
        stats['num_nconstraints'] = self.model.nconstraints ()

        if verbose:
            print ( "Considering %d activity windows" % (stats['num_acts']))
            print ( "Considering %d observation windows" % (stats['num_obs']))
            print ( "Considering %d dlink windows" % (stats['num_dlnks']))
            print ( "Considering %d xlink windows" % (stats['num_xlnks']))
            print ( 'self.model.nvariables ()')
            print ( self.model.nvariables ())
            print ( 'self.model.nobjectives ()')
            print ( self.model.nobjectives ())
            print ( 'self.model.nconstraints ()')
            print ( self.model.nconstraints ())

        return stats

    # def get_act_obs_subscripts(obs_winds,dlnk_winds_flat,xlnk_winds):
    #     lnk_obs_subscripts = []

    #     for sat_acts in 



    def filter_windows(self,obs_winds,dlnk_winds_flat,xlnk_winds_flat,num_sats):

        obs_winds_filtered = [[] for sat_indx in  range (num_sats)]
        dlink_winds_flat_filtered = [[] for sat_indx in  range (num_sats)]
        xlink_winds_flat_filtered = [[] for sat_indx in  range (num_sats)]

        for sat_indx in  range (num_sats):
            for wind in obs_winds[sat_indx]:
                if  wind.start >= self.planning_start_dt  and  wind.end  <= self.planning_end_obs_xlnk_dt:
                    obs_winds_filtered[sat_indx]. append ( wind)

            for wind in dlnk_winds_flat[sat_indx]:
                if  wind.start >= self.planning_start_dt  and  wind.end  <= self.planning_end_dlnk_dt:
                    dlink_winds_flat_filtered[sat_indx]. append ( wind)

            for wind in xlnk_winds_flat[sat_indx]:
                if  wind.start >= self.planning_start_dt  and  wind.end  <= self.planning_end_obs_xlnk_dt:
                    xlink_winds_flat_filtered[sat_indx]. append ( wind)

        return obs_winds_filtered, dlink_winds_flat_filtered, xlink_winds_flat_filtered

    def get_activity_structs(self,obs_winds_filt,dlnk_winds_filt,xlnk_winds_filt):
        """Get the structures representing activity subscripts for variables, parameters that are used in creating the pyomo model"""

        #  all activities are uniquely indexed. these structures keep track of those, and the mapping to activity objects
        all_acts_windids = set()
        all_obs_windids = set()
        all_dlnk_windids = set()
        all_xlnk_windids = set()
        all_act_objs_by_windid = {}
        capacity_by_act_windid = {}

        # all acts for a given sat
        sats_acts_set = [set() for sat_indx in range (self.num_sats)]

        # this is T^s - defines a unique time system for every satellite based on where the activities actually take place, as opposed to a finely granularized, monolithic list of times that could impose MILP constraints even when there are no activities happening during most of the times.
        all_link_times_by_sat_indx = {sat_indx: [] for sat_indx in range(self.num_sats)}  
        # this is used for producing T^s
        all_link_times_link_objs_by_sat_indx = {sat_indx: [] for sat_indx in range(self.num_sats)}  
        # this is a list of all the (s,o,t) in v_(s,o,t) in the formulation (where t is an index, starting at 0 for each sat)
        sat_obs_time_dv_subscripts = []
        # this is the actual time corresponding to the subscript stored in sat_obs_time_dv_subscripts
        time_by_sat_obs_time_dv_subscript = {}
        # keeps track of the v_(s,o,t) subscripts corresponding to "the start of" a given activity 
        sat_obs_time_dv_subscripts_by_link_obj = {}

        # making these sets because xlnks are repeated across xlnk.sat_indx, xlnk.xsat_indx, and want to use same data structure for both dlnks and xlnks
        dlnk_obs_windids = set()
        xlnk_obs_windids = set()

        obs_windids_by_dlnk_windid = {}
        dlnk_windids_by_obs_windid = {}
        obs_windids_by_xlnk_windid = {}
        xlnk_windids_by_obs_windid = {}

        # look for every dlnk that follows an obs. This is a potential window for downlinking data from the obs
        for sat_indx in range(self.num_sats):
            for dlnk in dlnk_winds_filt[sat_indx]:
                # grab the windid/obj mapping while we're here
                dlnk_windid = dlnk.window_ID
                all_act_objs_by_windid[dlnk_windid] = dlnk
                all_acts_windids.add(dlnk_windid)
                all_dlnk_windids.add(dlnk_windid)
                capacity_by_act_windid[dlnk_windid] = dlnk.data_vol
                obs_windids_by_dlnk_windid.setdefault(dlnk_windid, set())
                all_link_times_link_objs_by_sat_indx[sat_indx].append((dlnk.center,[dlnk]))
                sats_acts_set[sat_indx].add(dlnk)
                sat_obs_time_dv_subscripts_by_link_obj.setdefault(dlnk, [])

                for obs_sat in obs_winds_filt:
                    for obs in obs_sat:
                        # we'll do this obs storing again for every subsequent downlink, of course...it's fine, we'll survive
                        obs_windid = obs.window_ID
                        all_act_objs_by_windid[obs_windid] = obs
                        all_acts_windids.add(obs_windid)
                        all_obs_windids.add(obs_windid)
                        capacity_by_act_windid[obs_windid] = obs.data_vol
                        sats_acts_set[obs.sat_indx].add(obs)
                        
                        dlnk_windids_by_obs_windid.setdefault(obs_windid, set())

                        if dlnk.center > obs.center:
                            # make tuple of windids
                            dlnk_obs_windids.add((dlnk.window_ID,obs.window_ID))
                            obs_windids_by_dlnk_windid[dlnk_windid].add(obs_windid)
                            dlnk_windids_by_obs_windid[obs_windid].add(dlnk_windid)

        # look for every xlnk that follows an obs. This is a potential window for moving data from the obs
        # note: assuming we want to look at both xlnk directionalities here
        for sat_indx in range(self.num_sats):
            for xlnk in xlnk_winds_filt[sat_indx]:
                xlnk_windid = xlnk.window_ID
                all_act_objs_by_windid[xlnk_windid] = xlnk
                all_acts_windids.add(xlnk_windid)
                all_xlnk_windids.add(xlnk_windid)
                capacity_by_act_windid[xlnk_windid] = xlnk.data_vol
                # need set() here because xlnk windows are seen by sat_indcs on both sides of the link
                obs_windids_by_xlnk_windid.setdefault(xlnk_windid, set())
                # note that we can't use xlnk.sat_indx below below every xlnk wind has a sat_indx and an xsat_indx (sat_indx < xsat_indx) and it's incorrect to use sat_indx straight
                all_link_times_link_objs_by_sat_indx[sat_indx].append((xlnk.center,[xlnk]))
                sats_acts_set[sat_indx].add(xlnk)
                sat_obs_time_dv_subscripts_by_link_obj.setdefault(xlnk, [])

                for obs_sat in obs_winds_filt:
                    for obs in obs_sat:

                        obs_windid = obs.window_ID
                        xlnk_windids_by_obs_windid.setdefault(obs_windid, set())

                        if xlnk.center > obs.center:
                            # make tuple of windids
                            xlnk_obs_windids.add((xlnk.window_ID,obs.window_ID))
                            obs_windids_by_xlnk_windid[xlnk_windid].add(obs.window_ID)
                            xlnk_windids_by_obs_windid[obs_windid].add(xlnk_windid)

        # this update() call is fine because dlnk_windids, xlnk_windids are mutex
        obs_windids_by_lnk_windid = copy(obs_windids_by_dlnk_windid)
        obs_windids_by_lnk_windid.update(obs_windids_by_xlnk_windid)

        #  make sure they're sorted
        for time_link_objs_list in all_link_times_link_objs_by_sat_indx.values():
            # sort by time
            time_link_objs_list.sort(key=lambda x: x[0])


        all_link_times_link_objs_by_sat_indx_merged = {sat_indx: [] for sat_indx in range(self.num_sats)}  

        # deal with the fact that some activities have pretty much the same center time, and we'd like to merge those into a single timepoint in T_s rather than multiple, to avoid that extra hit in computational complexity
        for sat_indx,time_link_objs_list in all_link_times_link_objs_by_sat_indx.items():
            # grab the first time in the list
            # curr_time = time_link_objs_list[0][0]
            curr_time = time_link_objs_list[0][0]
            curr_link_objs = []

            for time_obj in time_link_objs_list:
                new_time = time_obj[0]
                new_link_objs = time_obj[1] # this should actually be only one link obj

                # if the new time is effectively the same as the old time...
                if new_time - curr_time < self.T_s_time_epsilon_td:
                    curr_link_objs += new_link_objs
                else:
                
                    all_link_times_link_objs_by_sat_indx_merged[sat_indx].append((curr_time,curr_link_objs))
                    all_link_times_by_sat_indx.setdefault(sat_indx,[])
                    all_link_times_by_sat_indx[sat_indx].append(curr_time)

                    curr_time = new_time
                    curr_link_objs = copy(new_link_objs) # copy the list just cuz it feels safer

            # flush whatever's left at the end
            all_link_times_link_objs_by_sat_indx_merged[sat_indx].append((curr_time,curr_link_objs))
            all_link_times_by_sat_indx.setdefault(sat_indx,[])
            all_link_times_by_sat_indx[sat_indx].append(curr_time)

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        # print("time_epsilon_tiny_diff_count")
        # print(time_epsilon_tiny_diff_count)

        # generate list of subscripts for v_(s,o,t), as well as the correspondence between activity windows and their entries in v_(s,o,t)
        for obs_sat in obs_winds_filt:
            for obs in obs_sat:
                for sat_indx in range(self.num_sats):
                    for time_index,(time,link_objs) in enumerate(all_link_times_link_objs_by_sat_indx_merged[sat_indx]):

                        if time > obs.center:
                            sat_obs_time_subscr = (sat_indx,obs.window_ID,time_index)
                            sat_obs_time_dv_subscripts.append(sat_obs_time_subscr)
                            time_by_sat_obs_time_dv_subscript[sat_obs_time_subscr] = time
                            for link_obj in link_objs:
                                # if sat_indx ==4 and obs.window_ID == 41 and link_obj.window_ID == 5303:
                                #     from circinus_tools import debug_tools
                                #     debug_tools.debug_breakpt()

                                sat_obs_time_dv_subscripts_by_link_obj[link_obj].append(sat_obs_time_subscr)


        sats_acts = []
        for sat_act_set in sats_acts_set:
            sat_act_list = list(sat_act_set)
            sat_act_list.sort(key= lambda a: a.center)
            sats_acts.append(sat_act_list)



        return sats_acts,all_acts_windids,all_obs_windids,all_dlnk_windids,all_xlnk_windids,all_act_objs_by_windid,capacity_by_act_windid,dlnk_obs_windids,xlnk_obs_windids,obs_windids_by_lnk_windid,obs_windids_by_dlnk_windid,obs_windids_by_xlnk_windid,dlnk_windids_by_obs_windid,xlnk_windids_by_obs_windid,all_link_times_by_sat_indx,sat_obs_time_dv_subscripts,time_by_sat_obs_time_dv_subscript,sat_obs_time_dv_subscripts_by_link_obj


    def make_model ( self,obs_winds,dlnk_winds_flat,xlnk_winds_flat, ecl_winds, verbose = False):
        # note this model only works for non-symmetric crosslink windows!

        model = pe.ConcreteModel()

        # filter the routes to make sure that  none of their activities fall outside the scheduling window
        obs_winds_filt,dlnk_winds_filt,xlnk_winds_filt = self.filter_windows(obs_winds,dlnk_winds_flat,xlnk_winds_flat,self.num_sats)

        self.ecl_winds = ecl_winds

        #  really useful code below!!!
        if verbose:
            pass

        ##############################
        #  Make indices/ subscripts
        ##############################

        try:
            # (sat_acts,
            #     sat_dlnks,
            #     all_acts_indcs,
            #     dmr_indcs_by_act,
            #     dv_by_act,
            #     all_acts_by_obj,
            #     all_acts_by_indx,
            #     obs_act_indcs,
            #     dmr_indcs_by_obs_act,
            #     dv_by_obs_act,
            #     link_act_indcs,
            #     dmr_indcs_by_link_act,
            #     dv_by_link_act) =  self.get_activity_structs(routes_flat)

            (sats_acts,
                all_acts_windids,
                all_obs_windids,
                all_dlnk_windids,
                all_xlnk_windids,
                all_act_objs_by_windid,
                capacity_by_act_windid,
                dlnk_obs_windids,
                xlnk_obs_windids,
                obs_windids_by_lnk_windid,
                obs_windids_by_dlnk_windid,
                obs_windids_by_xlnk_windid,
                dlnk_windids_by_obs_windid,
                xlnk_windids_by_obs_windid,
                all_link_times_by_sat_indx,
                sat_obs_time_dv_subscripts,
                time_by_sat_obs_time_dv_subscript,
                sat_obs_time_dv_subscripts_by_link_obj
                 ) = self.get_activity_structs(
                        obs_winds_filt,
                        dlnk_winds_filt,
                        xlnk_winds_filt
                    )

            self.all_acts_windids = all_acts_windids
            self.all_obs_windids = all_obs_windids
            self.all_dlnk_windids = all_acts_windids
            self.all_xlnk_windids = all_acts_windids
            self.all_act_objs_by_windid = all_act_objs_by_windid
            self.time_by_sat_obs_time_dv_subscript = time_by_sat_obs_time_dv_subscript
            # combine xlnk and dlnks windids by obs windid together
            self.lnk_windids_by_obs_windid = copy(xlnk_windids_by_obs_windid)
            for obs_windid,dlnk_windids in dlnk_windids_by_obs_windid.items():
                self.lnk_windids_by_obs_windid.setdefault(obs_windid, set())
                self.lnk_windids_by_obs_windid[obs_windid] = self.lnk_windids_by_obs_windid[obs_windid].union(dlnk_windids)


            # dmr_latency_sf_by_dmr_indx =  self.get_dmr_latency_score_factors(
            #     routes_flat,
            #     dmr_indcs_by_obs_act,
            #     self.latency_params
            # )

            # construct a set of dance cards for every satellite, 
            # each of which keeps track of all of the activities of satellite 
            # can possibly execute at any given time slice delta T. 
            # this is for constructing energy storage constraints
            # using resource_delta_t_s because this dancecard is solely for use in constructing resource constraints
            # note that these dancecards will baloon in size pretty quickly as the planning_end_dt increases. However most of the complexity they introduce is before planning_end_obs_xlnk_dt because that's the horizon where obs,xlnk actitivities are included. After that there should only be sparse downlinks
        #     es_act_dancecards = [Dancecard(self.planning_start_dt,self.planning_end_dt,self.resource_delta_t_s,item_init=None,mode='timestep') for sat_indx in range (self.num_sats)]
            
        #     for sat_indx in range (self.num_sats): 
        #         es_act_dancecards[sat_indx].add_winds_to_dancecard(sat_acts[sat_indx])
        #         es_act_dancecards[sat_indx].add_winds_to_dancecard(ecl_winds[sat_indx])

        #     # this is for data storage
        #     # for each sat/timepoint, we store a list of all those data multi routes that are storing data on the sat at that timepoint
        #     ds_route_dancecards = [Dancecard(self.planning_start_dt,self.planning_end_dt,self.resource_delta_t_s,item_init=None,mode='timepoint') for sat_indx in range (self.num_sats)]
            
        #     # add data routes to the dancecard
        #     for dmr_indx,dmr in enumerate(routes_flat):
        #         # list of type routing_objects.SatStorageInterval
        #         dmr_ds_intervs = dmr.get_data_storage_intervals()

        #         for interv in dmr_ds_intervs:
        #             # store the dmr object at this timepoint
        #             ds_route_dancecards[interv.sat_indx].add_item_in_interval(dmr_indx,interv.start,interv.end)

        except IndexError:
            raise RuntimeWarning('sat_indx out of range. Are you sure all of your input files are consistent? (including pickles)')        
        # self.dmr_indcs_by_act = dmr_indcs_by_act
        # self.all_acts_indcs = all_acts_indcs
        # self.obs_act_indcs = obs_act_indcs
        # self.link_act_indcs = link_act_indcs
        # self.all_acts_by_obj = all_acts_by_obj
        # self.all_acts_by_indx = all_acts_by_indx

        # # these subscripts probably should've been done using the unique IDs for the objects, rather than their arbitrary locations within a list. Live and learn, hélas...

        # #  subscript for each dmr (data multi route) p  (p index is a hold-over from when I called them paths)
        # model.dmrs = pe.Set(initialize= range(len(routes_flat)))
        
        #  subscript for each activity a
        model.act_windids = pe.Set(initialize= all_acts_windids)
        #  subscript for each obs o
        model.obs_windids = pe.Set(initialize= all_obs_windids)
        #  subscript for each link
        model.lnk_windids = pe.Set(initialize= all_dlnk_windids.union(all_xlnk_windids))
        # every act-obs within that act combination
        model.lnk_obs_subscripts = pe.Set(initialize= dlnk_obs_windids.union(xlnk_obs_windids),dimen=2)
        #  subscript for each satellite
        model.sat_indcs = pe.Set(initialize= range(self.num_sats))
        #  subscript for each s,o,t in v_(s,o,t)
        # t in model.sat_obs_time_dv_subscripts is index into all_link_times_by_sat_indx
        model.sat_obs_time_dv_subscripts = pe.Set(initialize= sat_obs_time_dv_subscripts,dimen=3)


        # # timepoints is the indices, which starts at 0 
        # #  NOTE: we assume the same time system for every satellite
        # self.es_time_getter_dc = es_act_dancecards[0]
        # es_num_timepoints = es_act_dancecards[0].num_timepoints
        # model.es_timepoint_indcs = pe.Set(initialize=  self.es_time_getter_dc.get_tp_indcs())

        # self.ds_time_getter_dc = ds_route_dancecards[0]
        # model.ds_timepoint_indcs = pe.Set(initialize=  self.ds_time_getter_dc.get_tp_indcs())

        # #  unique indices for observation and link acts
        # model.obs_acts = pe.Set(initialize= obs_act_indcs)
        # model.link_acts = pe.Set(initialize= link_act_indcs)

        # if self.solver_name == 'gurobi' or self.solver_name == 'cplex':
        #     int_feas_tol = self.solver_params['integer_feasibility_tolerance']
        # else:
        #     raise NotImplementedError

        # for p,lat_sf in dmr_latency_sf_by_dmr_indx.items():        
        #     if lat_sf > int_feas_tol*self.big_M_lat:
        #         raise RuntimeWarning('big_M_lat (%f) is not large enough for latency score factor %f and integer feasibility tolerance %f (dmr index %d)'%(self.big_M_lat,lat_sf,int_feas_tol,p))

        # for act_obj in all_acts_by_obj.keys():
        #     if 2*(act_obj.end-act_obj.start).total_seconds() > self.big_M_act_t_dur_s:
        #         raise RuntimeWarning('big_M_act_t_dur_s (%f) is not large enough for act of duration %s and integer feasibility tolerance %f (act string %s)'%(self.big_M_act_t_dur_s,act_obj.end-act_obj.start,int_feas_tol,act_obj))
        #     # if 2*(act_obj.end-act_obj.start).total_seconds() > (1-int_feas_tol) *  self.big_M_act_t_dur_s:
        #         # raise RuntimeWarning('big_M_act_t_dur_s (%f) is not large enough for act of duration %s and integer feasibility tolerance %f (act string %s)'%(self.big_M_act_t_dur_s,act_obj.end-act_obj.start,int_feas_tol,act_obj))

        #     if 2*act_obj.data_vol > (1-int_feas_tol) * self.big_M_act_dv:
        #         raise RuntimeWarning('big_M_act_dv (%f) is not large enough for act of dv %f and integer feasibility tolerance %f (act string %s)'%(self.big_M_act_dv,act_obj.data_vol,int_feas_tol,act_obj))
                


        ##############################
        #  Make parameters
        ##############################

        # model.par_min_as_route_dv = pe.Param (initialize=self.min_as_route_dv)
        # model.par_obs_dv = pe.Param(model.obs_acts,initialize =dv_by_obs_act)
        model.par_total_obs_dv = sum(capacity_by_act_windid[o] for o in all_obs_windids)
        # model.par_link_dv = pe.Param(model.link_acts,initialize =dv_by_link_act)
        model.par_act_capacity = pe.Param(model.act_windids,initialize =capacity_by_act_windid)
        # #  data volume for each data multi-route
        # model.par_dmr_dv = pe.Param(model.dmrs,initialize ={ dmr_indx: dmr.data_vol for dmr_indx,dmr in enumerate (routes_flat)})
        # #  data volume for each activity in each data multi-route

        # # todo: will probably have to fix this to not be the full set multiplication of model.dmrs, model.acts - can use dmr_indcs_by_act
        # model.par_dmr_act_dv = pe.Param(
        #     model.dmrs,
        #     model.acts,
        #     initialize = { (dmr_indx,self.all_acts_by_obj[act]): 
        #         dmr.data_vol_for_wind(act) for dmr_indx,dmr in enumerate (routes_flat) for act in dmr.get_winds()
        #     }
        # )

        # # each of these is essentially a dictionary indexed by link or obs act indx, with  the value being a list of dmr indices that are included within that act
        # # these are valid indices into model.dmrs
        # model.par_dmr_subscrs_by_link_act = pe.Param(model.link_acts,initialize =dmr_indcs_by_link_act)
        # model.par_dmr_subscrs_by_obs_act = pe.Param(model.obs_acts,initialize =dmr_indcs_by_obs_act)
        # model.par_dmr_subscrs_by_act = pe.Param(model.acts,initialize =dmr_indcs_by_act)


        # if self.energy_unit == "Wh":
        #     model.par_resource_delta_t = pe.Param (initialize= self.resource_delta_t_s/3600)
        # else: 
        #     raise NotImplementedError
        # model.par_sats_estore_initial = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_init_estate_Wh)})
        # model.par_sats_estore_min = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_emin_Wh)})
        # model.par_sats_estore_max = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_emax_Wh)})
        # model.par_sats_edot_by_act = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_edot_by_act_W)})

        # model.par_sats_dstore_min = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_dmin_Mb)})
        # model.par_sats_dstore_max = pe.Param ( model.sats,initialize= { i: item for i,item in enumerate (self.sats_dmax_Mb)})

        ##############################
        #  Make variables
        ##############################

        # v_a, variables [1]
        model.var_act_dv_utilization  = pe.Var (model.act_windids, within = pe.NonNegativeReals)
        # v_(a,o), variables [2]
        model.var_lnk_obs_dv_utilization  = pe.Var (model.lnk_obs_subscripts,  within = pe.NonNegativeReals)
        # v_(s,o,t), variables [3]
        model.var_sat_obs_time_dv  = pe.Var (model.sat_obs_time_dv_subscripts,  within = pe.NonNegativeReals)




        # # activity utilization variable indicating how much of an activity's capacity is used [1]
        # model.var_activity_utilization  = pe.Var (model.acts, bounds =(0,1))
        # # dmr utilization variable indicating how much of a dmr's capacity is used [2]
        # model.var_dmr_utilization  = pe.Var (model.dmrs, bounds =(0,1))
        # #  indicator variables for whether or not dmrs [3] and activities [4] have been chosen
        # model.var_dmr_indic  = pe.Var (model.dmrs, within = pe.Binary)
        # # model.var_dmr_indic  = pe.Var (model.dmrs, bounds =(0,1))
        # model.var_act_indic  = pe.Var (model.acts, within = pe.Binary)
        # # model.var_act_indic  = pe.Var (model.acts, bounds =(0,1))

        
        # # satellite energy storage
        # model.var_sats_estore  = pe.Var (model.sats,  model.es_timepoint_indcs,  within = pe.NonNegativeReals)

        # # satellite data storage (data buffers)
        # model.var_sats_dstore  = pe.Var (model.sats,  model.ds_timepoint_indcs,  within = pe.NonNegativeReals)

        # model.var_dmr_latency_sf_by_obs_indx = pe.Var (model.obs_acts,  within = pe.NonNegativeReals)
        
        # allow_act_timing_constr_violations = False
        # if allow_act_timing_constr_violations:
        #     print('allow_act_timing_constr_violations is True')

        # #  variables for handling the allowance of inter-activity timing constraint violations. these are only generated if allow_act_timing_constr_violations is True
        # model.var_intra_sat_act_constr_violations = pe.VarList()
        # model.var_inter_sat_act_constr_violations = pe.VarList()

        # #  stores all of the lower bounds of the constraint violation variables, for use in normalization for objective function
        # min_var_intra_sat_act_constr_violation_list = [] 
        # min_var_inter_sat_act_constr_violation_list = [] 

        ##############################
        #  Make constraints
        ##############################

        # TODO: renumber  these with the final numbering

        # note that the observations show up within model.acts as well, so we also constraint route scheduled DV by the real available DV from each observation
        def c1_rule( model,a):
            return ( model.var_act_dv_utilization[a] <= model.par_act_capacity[a])
        model.c1 =pe.Constraint ( model.act_windids,  rule=c1_rule)

        def c1b_rule( model,a):
            return model.var_act_dv_utilization[a] >= sum(model.var_lnk_obs_dv_utilization[a,o] for o in obs_windids_by_lnk_windid[a])
        model.c1b =pe.Constraint ( model.lnk_windids,  rule=c1b_rule)

        # c2
        model.c2  = pe.ConstraintList()
        for a in model.lnk_windids:
            for o in obs_windids_by_lnk_windid[a]:
                model.c2.add(model.var_lnk_obs_dv_utilization[a,o] <= model.var_act_dv_utilization[o])


        def c4_rule( model,o):
            return model.var_act_dv_utilization[o] == sum(model.var_lnk_obs_dv_utilization[d,o] for d in dlnk_windids_by_obs_windid[o])
        model.c4 =pe.Constraint ( model.obs_windids,  rule=c4_rule)

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        # c5
        model.c5  = pe.ConstraintList()
        for o in model.obs_windids:
            for s in model.sat_indcs:
                t_o = all_act_objs_by_windid[o].center

                for t_indx,t in enumerate(all_link_times_by_sat_indx[s]):
                    if not t > t_o:
                        continue

                    # sum of all the dv terms up to time t
                    dv_sum = 0

                    # if this is the observing sat, account for the data given to the sat by the obs itself
                    if s == all_act_objs_by_windid[o].sat_indx:
                        dv_sum += model.var_act_dv_utilization[o]

                    # all of the acts that occur on s after obs (t_o) and before t
                    # note the inequalities here. First inequality - don't count the obs. Second: want v_(s,o,t) to be constrained by what happened previously. If this were an inequality, would say "the dv available at the time of an xlnk/dlnk is limited by dv contributions of all previous xlnks and the dv contribution of THIS XLNK. That's unbounded, could generate infinite dv"
                    sat_acts_slice = [act for act in sats_acts[s] if (
                        (type(act) == DlnkWindow or type(act) == XlnkWindow) and
                        (act.center > t_o and act.center < t)
                    )]
                    
                    for lnk in sat_acts_slice:
                        if type(lnk) == DlnkWindow:
                            # we lose dv in a dlnk window

                            dv_sum -= model.var_lnk_obs_dv_utilization[lnk.window_ID,o]
                        elif type(lnk) == XlnkWindow:
                            # this code won't work if we're assuming symmetric xlnk windows
                            assert(not lnk.symmetric)

                            # if receiving, we're adding dv in that xlnk window
                            if lnk.is_rx(s):
                                dv_sum += model.var_lnk_obs_dv_utilization[lnk.window_ID,o]
                            # if transmitting, we're losing dv in that xlnk window
                            else:
                                dv_sum -= model.var_lnk_obs_dv_utilization[lnk.window_ID,o]

                    # note using inequality to reduce computation time (maybe)
                    model.c5.add(model.var_sat_obs_time_dv[s,o,t_indx] == dv_sum)

                    # if s==4 and o == 41 and t_indx == 205:
                    # if s==4 and o == 64:
                    #     from circinus_tools import debug_tools
                    #     debug_tools.debug_breakpt()    

        # c6
        model.c6  = pe.ConstraintList()
        for s in model.sat_indcs:
            for lnk in sats_acts[s]:
                # verify lnk is actually a link act
                if not (type(lnk) == DlnkWindow or type(lnk) == XlnkWindow):
                    continue

                # this loop handles all observations that a link might be able to route
                for sat_obs_time_dv_subscript in sat_obs_time_dv_subscripts_by_link_obj[lnk]:
                    
                    # if receiving, we're not constrained by any dv availability
                    # need to consider this here because the sat_obs_time_dv_subscripts_by_link_obj mapping holds both tx and rx directions for every link
                    sat_indx_tmp = sat_obs_time_dv_subscript[0]
                    if type(lnk) == XlnkWindow: 
                        # this code won't work if we're assuming symmetric xlnk windows
                        assert(not lnk.symmetric)
                        if lnk.is_rx(sat_indx_tmp):
                                continue

                    # quick sanity check that it's the right satellite while we're here
                    # assert(s == sat_obs_time_dv_subscript[0])
                    o = sat_obs_time_dv_subscript[1]
                    model.c6.add( model.var_lnk_obs_dv_utilization[lnk.window_ID,o] <= model.var_sat_obs_time_dv[sat_obs_time_dv_subscript])

                    # if lnk.window_ID == 6735:
                    #     print(model.c6[len(model.c6)].expr.to_string())
                    # from circinus_tools import debug_tools
                    # debug_tools.debug_breakpt()
                    


        # # keep track of this, so we know to warn user about default transition time usage
        # #  note: this is not a parameter! it's merely helpful for reporting warnings
        # used_default_transition_time = False

        # def gen_inter_act_constraint(var_list,constr_list,transition_time_req,act1,act2,act1_uindx,act2_uindx):

        #     var_constr_violation = None
        #     min_constr_violation = None

        #     center_time_diff = (act2.center - act1.center).total_seconds()
        #     time_adjust_1 = model.par_act_dv[act1_uindx]*model.var_activity_utilization[act1_uindx]/2/act1.ave_data_rate
        #     time_adjust_2 = model.par_act_dv[act2_uindx]*model.var_activity_utilization[act2_uindx]/2/act2.ave_data_rate
        #     max_time_adjust_1 = model.par_act_dv[act1_uindx]*1/2/act1.ave_data_rate
        #     max_time_adjust_2 = model.par_act_dv[act2_uindx]*1/2/act2.ave_data_rate

        #     if not allow_act_timing_constr_violations:
        #         #  if the activities overlap in center time (including  transition time), then it's not possible to have sufficient transition time between them.  only allow one
        #         if (act2.center - act1.center).total_seconds() <= transition_time_req:
        #             constr = model.var_act_indic[act1_uindx]+ model.var_act_indic[act2_uindx] <= 1

        #         # If they don't overlap in center time, but they do overlap to some amount, then we need to constrain their end and start times to be consistent with one another
        #         else:
        #             M = max(act1.duration,act2.duration).total_seconds()
        #             constr_disable_1 = M*(1-model.var_act_indic[act1_uindx])
        #             constr_disable_2 = M*(1-model.var_act_indic[act2_uindx])
        
        #             constr = center_time_diff - time_adjust_1 - time_adjust_2 + constr_disable_1 + constr_disable_2 >= transition_time_req

        #     else:
        #         # time_adjust_N can go as low as zero, so the constraint violation can be this at its lowest
        #         min_constr_violation = center_time_diff - max_time_adjust_1 - max_time_adjust_2 - transition_time_req

        #         assert(min_constr_violation < 0)

        #         # deal with adding a variable to represent this constraint violation. ( I hate that I have to do it this way, deep within this function, but it seems like the approach you have to use for dynamically generating variable lists in pyomo, bizarrely. refer to: https://projects.coin-or.org/Coopr/browser/pyomo/trunk/pyomo/core/base/var.py?rev=11067, https://groups.google.com/forum/#!topic/pyomo-forum/modS1VkPxW0
        #         var_list.add()
        #         var_constr_violation = var_list[len(var_list)]

        #         #  bounds on range of constraint violation variable

        #         # bound minimum of the constraint violation (where both activities have 1.0 utilization), to keep the problem formulation as tight as possible
        #         var_constr_violation.setlb(min_constr_violation)
        #         # want the constraint violation only to go to zero at its maximum, because we don't want to reward times where there is no constraint violation, only penalize
        #         var_constr_violation.setub(0)

        #         #  the actual time constraint that bounds the constraint violation
        #         constr = center_time_diff - time_adjust_1 - time_adjust_2 - transition_time_req >= var_constr_violation

        #     return constr, var_constr_violation, min_constr_violation


        # #  intra-satellite activity overlap constraints [4],[5],[5b]
        # #  well, 5B is activity minimum time duration
        # # model.c4  = pe.ConstraintList()  # c5 now holds c4 constraints
        # model.c5  = pe.ConstraintList() # this now contains all of the activity overlap constraints
        # model.c5b  = pe.ConstraintList()
        # model.intra_sat_act_constr_bounds  = pe.ConstraintList()
        # self.intra_sat_act_constr_violation_acts_list = []
        # for sat_indx in range (self.num_sats):
        #     num_sat_acts = len(sat_acts[sat_indx])
        #     for  first_act_indx in  range (num_sat_acts):

        #         act1 = sat_acts[sat_indx][first_act_indx]
        #         # get the unique index into model.acts
        #         act1_uindx = all_acts_by_obj[act1]
        #         length_1 = model.par_act_dv[act1_uindx]*model.var_activity_utilization[act1_uindx]/act1.ave_data_rate
        #         model.c5b.add( length_1 >= model.var_act_indic[act1_uindx] * self.min_act_duration_s[type(act1)])
                
        #         for  second_act_indx in  range (first_act_indx+1,num_sat_acts):
        #             act2 = sat_acts[sat_indx][ second_act_indx]
        #             # get the unique index into model.acts
        #             act2_uindx = all_acts_by_obj[act2]

        #             # act list should be sorted
        #             assert(act2.center >= act1.center)

        #             # get the transition time requirement between these activities
        #             try:
        #                 transition_time_req = self.act_transition_time_map[("intra-sat",type(act1),type(act2))]
        #             # if not explicitly specified, go with default transition time requirement
        #             except KeyError:
        #                 used_default_transition_time = True
        #                 transition_time_req = self.act_transition_time_map["default"]

        #             # if there is enough transition time between the two activities, no constraint needs to be added
        #             #  note that we are okay even if for some reason Act 2 starts before Act 1 ends, because time deltas return negative total seconds as well
        #             if (act2.start - act1.end).total_seconds() >= transition_time_req:
        #                 #  don't need to do anything,  continue on to next activity pair
        #                 continue

        #             else:

        #                 constr, var_constr_violation, min_constr_violation = gen_inter_act_constraint(
        #                     model.var_intra_sat_act_constr_violations,
        #                     model.intra_sat_act_constr_bounds,
        #                     transition_time_req,
        #                     act1,
        #                     act2,
        #                     act1_uindx,
        #                     act2_uindx
        #                 )

        #                 #  add the constraint, regardless of whether or not it's a "big M" constraint, or a constraint violation constraint - they're handled the same
        #                 model.c5.add( constr )

        #                 #  if it's a constraint violation constraint, then we have a variable to deal with
        #                 if not min_constr_violation is None:
        #                     min_var_intra_sat_act_constr_violation_list.append(min_constr_violation)
        #                     self.intra_sat_act_constr_violation_acts_list.append((act1,act2))


        # # inter-satellite downlink overlap constraints [9],[10]
        # # model.c9  = pe.ConstraintList() # c10 now holds c9 constraints
        # model.c10  = pe.ConstraintList()
        # model.inter_sat_act_constr_bounds  = pe.ConstraintList()
        # self.inter_sat_act_constr_violation_acts_list = []
        # for sat_indx in range (self.num_sats):
        #     num_sat_acts = len(sat_dlnks[sat_indx])
            
        #     for other_sat_indx in range (self.num_sats):
        #         if other_sat_indx == sat_indx:
        #             continue

        #         num_other_sat_acts = len(sat_dlnks[other_sat_indx])

        #         for  sat_act_indx in  range (num_sat_acts):

        #             act1 = sat_dlnks[sat_indx][sat_act_indx]
        #             # get the unique index into model.acts
        #             act1_uindx = all_acts_by_obj[act1]
                    
        #             for  other_sat_act_indx in  range (num_other_sat_acts):
        #                 act2 = sat_dlnks[other_sat_indx][other_sat_act_indx]
        #                 # get the unique index into model.acts
        #                 act2_uindx = all_acts_by_obj[act2]

        #                 assert(type(act1) == DlnkWindow and type(act2) == DlnkWindow)

        #                 # this line is pretty important - only consider overlap if they're looking at the same GS. I forgot to add this before and spent days wondering why the optimization process was progressing so slowly (hint: it's really freaking constrained and there's not much guidance for finding a good objective value if no downlink can overlap in time with any other downlink)
        #                 if act1.gs_indx != act2.gs_indx:
        #                     continue

        #                 # we're considering windows across satellites, so they could be out of order temporally. These constraints are only valid if act2 is after act1 (center time). Don't worry, as we loop through satellites, we consider both directions (i.e. act1 and act2 will be swapped in another iteration, and we'll get past this check and impose the required constraints)
        #                 if (act2.center - act1.center).total_seconds() < 0:
        #                     continue


        #                 # get the transition time requirement between these activities
        #                 try:
        #                     transition_time_req = self.act_transition_time_map[("inter-sat",DlnkWindow,DlnkWindow)]
        #                 # if not explicitly specified, go with default transition time requirement
        #                 except KeyError:
        #                     used_default_transition_time = True
        #                     transition_time_req = self.act_transition_time_map["default"]

        #                 # if there is enough transition time between the two activities, no constraint needs to be added
        #                 #  note that we are okay even if for some reason Act 2 starts before Act 1 ends, because time deltas return negative total seconds as well
        #                 if (act2.start - act1.end).total_seconds() >= transition_time_req:
        #                     #  don't need to do anything,  continue on to next activity pair
        #                     continue

        #                 else:
        #                     constr, var_constr_violation, min_constr_violation = gen_inter_act_constraint(
        #                         model.var_inter_sat_act_constr_violations,
        #                         model.inter_sat_act_constr_bounds,
        #                         transition_time_req,
        #                         act1,
        #                         act2,
        #                         act1_uindx,
        #                         act2_uindx
        #                     )

        #                     #  add the constraint, regardless of whether or not it's a "big M" constraint, or a constraint violation constraint - they're handled the same
        #                     model.c10.add( constr )

        #                     #  if it's a constraint violation constraint, then we have a variable to deal with
        #                     if not min_constr_violation is None:
        #                         # model.var_inter_sat_act_constr_violations.add(var_constr_violation)
        #                         min_var_inter_sat_act_constr_violation_list.append(min_constr_violation)
        #                         self.inter_sat_act_constr_violation_acts_list.append((act1,act2))


        # if verbose:
        #     if used_default_transition_time:
        #         print('\nWarning: used default transition time for inter- or intra- satellite activity timing constraints\n')


        # #  energy constraints [6]
        # model.c6  = pe.ConstraintList()
        # for sat_indx in range (self.num_sats): 

        #     # tp_indx serves as an index into the satellite activity dance cards
        #     for tp_indx in model.es_timepoint_indcs:
        #         #  constraining first time step to initial energy storage
        #         #  continue for loop afterwards because no geq/leq constraints needed for this index
        #         if tp_indx == 0:
        #             model.c6.add( model.var_sats_estore[sat_indx,tp_indx] ==  model.par_sats_estore_initial[sat_indx])
        #             continue 

        #         #  minimum and maximum storage constraints
        #         model.c6.add( model.var_sats_estore[sat_indx,tp_indx] >= model.par_sats_estore_min[sat_indx])
        #         model.c6.add( model.var_sats_estore[sat_indx,tp_indx] <= model.par_sats_estore_max[sat_indx])

        #         if self.enforce_energy_storage_constr:
        #             # determine activity energy consumption
        #             charging = True
        #             activity_delta_e = 0 
        #             #  get the activities that were active during the time step immediately preceding time point
        #             activities = es_act_dancecards[sat_indx].get_objects_at_ts_pre_tp_indx(tp_indx)
        #             # activities may be none if nothing is happening at timestep, to minimize RAM usage
        #             if activities:
        #                 for act in activities:
        #                     #  if this is a "standard activity" that we can choose to perform or not
        #                     if type(act) in self.standard_activities:
        #                         act_uindx = all_acts_by_obj[act]
        #                         activity_delta_e += (
        #                             model.par_sats_edot_by_act[sat_indx][self.act_type_map[type(act)]] 
        #                             * model.var_activity_utilization[act_uindx]
        #                             * model.par_resource_delta_t
        #                         )

        #                     if type(act) == XlnkWindow:
        #                         act_uindx = all_acts_by_obj[act]

        #                         if self.use_symmetric_xlnk_windows:
        #                             xlnk_edot = model.par_sats_edot_by_act[sat_indx]['xlnk-tx']
        #                         # note that in the case where we're not using symmetric xlnk windows, both satellites have a copy of each unidirectional window in their list of activity windows, so it'll be added appropriately at each unique sat_indx
        #                         else:
        #                             is_rx = act.is_rx(sat_indx)
        #                             if is_rx: 
        #                                 xlnk_edot = model.par_sats_edot_by_act[sat_indx]['xlnk-rx']
        #                             else: 
        #                                 xlnk_edot = model.par_sats_edot_by_act[sat_indx]['xlnk-tx']

        #                         activity_delta_e += (
        #                             xlnk_edot 
        #                             * model.var_activity_utilization[act_uindx]
        #                             * model.par_resource_delta_t
        #                         )

        #                     #  if the satellite is not in sunlight then we can't charge
        #                     elif type(act) == EclipseWindow:
        #                         charging = False

        #             # add in charging energy contribution ( if possible)
        #             charging_delta_e = model.par_sats_edot_by_act[sat_indx]['charging']*model.par_resource_delta_t if charging else 0

        #             #  base-level satellite energy usage (not including additional activities)
        #             base_delta_e = model.par_sats_edot_by_act[sat_indx]['base']*model.par_resource_delta_t

        #             # maximum bound of energy at current time step based on last time step
        #             model.c6.add( model.var_sats_estore[sat_indx,tp_indx] <= 
        #                 model.var_sats_estore[sat_indx,tp_indx-1]
        #                 + activity_delta_e
        #                 + charging_delta_e
        #                 + base_delta_e
        #             )

        #             # maximum bound of energy at current time step based on last time step
        #             model.c6.add( model.var_sats_estore[sat_indx,tp_indx] >= 
        #                 model.var_sats_estore[sat_indx,tp_indx-1]
        #                 + activity_delta_e
        #                 + base_delta_e
        #             )


        # #  data storage constraints [7]
        # model.c7  = pe.ConstraintList()
        # for sat_indx in range (self.num_sats): 

        #     # tp_indx serves as an index into the satellite data storage dance cards
        #     for tp_indx in model.ds_timepoint_indcs:

        #         # todo: add in an intital data volume value?
        #         #  maximum storage constraints
        #         model.c7.add( model.var_sats_dstore[sat_indx,tp_indx] <= model.par_sats_dstore_max[sat_indx])

        #         if self.enforce_data_storage_constr:
        #             routes_storing = ds_route_dancecards[sat_indx][tp_indx]

        #             # todo: this may be too slow to use == below. Change it back to >= and use a better approach to extract real data storage values in extract_resource_usage below? Leaving == for now because it means I can use these vars to extract output data usage values
        #             # constrain the minimum data storage at this tp by the amount of data volume being buffered by each route passing through the sat
        #             if routes_storing:
        #                 model.c7.add( model.var_sats_dstore[sat_indx,tp_indx] == sum(model.par_dmr_dv[p]*model.var_dmr_utilization[p] for p in routes_storing))
        #             else:
        #                 model.c7.add( model.var_sats_dstore[sat_indx,tp_indx] == 0)



        # #  observation latency score factor constraints [8]
        # model.c8  = pe.ConstraintList()
        # for obs_indx in model.obs_acts:
        #     dmrs_obs = model.par_dmr_subscrs_by_obs_act[obs_indx]

        #     #  sort the latency score factors for all the dmrs for this observation in increasing order -  important for constraint construction
        #     dmrs_obs.sort(key= lambda p: dmr_latency_sf_by_dmr_indx[p])

        #     num_dmrs_obs = len(dmrs_obs)
        #     #  initial constraint -  score factor for this observation will be equal to zero if no dmrs for this obs were chosen
        #     model.c8.add( model.var_dmr_latency_sf_by_obs_indx[obs_indx] <= 0 + self.big_M_lat * sum(model.var_dmr_indic[p] for p in dmrs_obs) )
            
        #     for dmr_obs_indx in range(num_dmrs_obs):
        #         #  add constraint that score factor for observation is less than or equal to the score factor for this dmr_obs_indx, plus any big M terms for any dmrs with larger score factors.
        #         #  what this does is effectively disable the constraint for the score factor for this dmr_obs_indx if any higher score factor dmrs were chosen 
        #         model.c8.add( model.var_dmr_latency_sf_by_obs_indx[obs_indx] <= 
        #             dmr_latency_sf_by_dmr_indx[dmrs_obs[dmr_obs_indx]] + 
        #             self.big_M_lat * sum(model.var_dmr_indic[p] for p in dmrs_obs[dmr_obs_indx+1:num_dmrs_obs]) )

        #         #  note: use model.c8[indx].expr.to_string()  to print out the constraint in a human readable form
        #         #                ^ USES BASE 1 INDEXING!!! WTF??
                

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        ##############################
        #  Make objective
        ##############################


        # #  determine which time points to use as "spot checks" on resource margin. These are the points that will be used in the objective function for maximizing resource margin
        # timepoint_spacing = ceil(es_num_timepoints/self.resource_margin_obj_num_timepoints)
        # # need to turn the generator into a list for slicing
        # #  note: have to get the generator again
        # decimated_tp_indcs = list(self.es_time_getter_dc.get_tp_indcs())[::timepoint_spacing]
        # rsrc_norm_f = len(decimated_tp_indcs) * len(model.sats)

        def obj_rule(model):
            total_dv_term = self.obj_weights['obs_dv'] * 1/model.par_total_obs_dv * sum(model.var_act_dv_utilization[o] for o in model.obs_windids) 

            # latency_term = self.obj_weights['route_latency'] * 1/len(model.obs_acts) * sum(model.var_dmr_latency_sf_by_obs_indx[o] for o in model.obs_acts)
            
            # energy_margin_term = self.obj_weights['energy_storage'] * 1/rsrc_norm_f * sum(model.var_sats_estore[sat_indx,tp_indx]/model.par_sats_estore_max[sat_indx] for tp_indx in decimated_tp_indcs for sat_indx in model.sats)

            # if len(min_var_inter_sat_act_constr_violation_list) > 0:
            #     inter_sat_act_constr_violations_term = self.obj_weights['inter_sat_act_constr_violations'] * 1/sum(min_var_inter_sat_act_constr_violation_list) * sum(model.var_inter_sat_act_constr_violations[indx] for indx in range(1,len(model.var_inter_sat_act_constr_violations)+1))
            # else:
            #     inter_sat_act_constr_violations_term = 0

            # if len(min_var_intra_sat_act_constr_violation_list) > 0:
            #     intra_sat_act_constr_violations_term = self.obj_weights['intra_sat_act_constr_violations'] * 1/sum(min_var_intra_sat_act_constr_violation_list) * sum(model.var_intra_sat_act_constr_violations[indx] for indx in range(1,len(model.var_inter_sat_act_constr_violations)+1))
            # else:
            #     intra_sat_act_constr_violations_term = 0

            # from circinus_tools import debug_tools
            # debug_tools.debug_breakpt()

            # return total_dv_term + latency_term + energy_margin_term - inter_sat_act_constr_violations_term - intra_sat_act_constr_violations_term
            return total_dv_term
            
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

    def print_sol(self):
        for v in self.model.component_objects(pe.Var, active=True):
            if str (v) =='var_act_dv_utilization': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    # index is actually a window id
                    act = self.all_act_objs_by_windid[index]
                    if type(act) == ObsWindow:
                        print (" %d %0.3f   %s"%(index, val, act))

            if str (v) =='var_lnk_obs_dv_utilization': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    # index is actually a window id
                    # if type(self.all_act_objs_by_windid[index[0]]) == DlnkWindow:
                    if val>0.001:
                        print (" %s %0.3f   %s %s"%(index, val, self.all_act_objs_by_windid[index[1]],self.all_act_objs_by_windid[index[0]]))


            if str (v) =='var_sat_obs_time_dv': 
                sorted_var_sat_obs_time_dv = [[] for sat_indx in range(self.num_sats)]

                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value

                    # index is (sat_indx,obs windid, time index)
                    sat_indx = index[0]
                    obs_windid = index[1]
                    time_indx = index[2]
                    sorted_var_sat_obs_time_dv[sat_indx].append(index)

                    # if sat_indx == 3 and  obs_windid == 41:
                    #     print (" %s %0.3f   %s"%(index, val, self.time_by_sat_obs_time_dv_subscript[index]))

                for sat_indx in range(self.num_sats):
                    # sort by time index
                    sorted_var_sat_obs_time_dv[sat_indx].sort(key= lambda subscr: subscr[2])

                    for index in sorted_var_sat_obs_time_dv[sat_indx]:
                        val  = varobject[index].value
                        obs_windid = index[1]
                        # if sat_indx == 4 and  obs_windid == 41:
                        if sat_indx == 4 and obs_windid == 64:
                            print (" %s %0.3f   %s"%(index, val, self.time_by_sat_obs_time_dv_subscript[index]))

    def fabricate_routes(self):

        found_routes = []
        dr_uid = 0
        for o in self.model.obs_windids:

            # if not enough dv from this obs, then don't construct routes for it
            if not (pe.value(self.model.var_act_dv_utilization[o]) > self.min_as_route_dv):
                continue

            remaining_obs_dv = pe.value(self.model.var_act_dv_utilization[o])

            obs = self.all_act_objs_by_windid[o]

            # construct objects representing all the links available for this obs, and the available capacity on each of them for routing data from this obs
            lnk_windids_o = self.lnk_windids_by_obs_windid[o]
            remaining_dv_o_by_lnk = {}
            # keeps track of links originating from a given sat index
            lnks_o_by_sat_indx = [[] for sat_indx in range(self.num_sats)]

            for a in lnk_windids_o:
                lnk_dv_o = pe.value(self.model.var_lnk_obs_dv_utilization[a,o])

                if lnk_dv_o > self.dv_epsilon:
                    lnk = self.all_act_objs_by_windid[a]
                    remaining_dv_o_by_lnk[lnk] = lnk_dv_o

                    if type(lnk) == DlnkWindow:
                        lnks_o_by_sat_indx[lnk.sat_indx].append(lnk)
                    elif type(lnk) == XlnkWindow:
                        lnks_o_by_sat_indx[lnk.tx_sat].append(lnk)
                    else:
                        # sanity check
                        raise RuntimeWarning("Should only find dlnks and xlnks here")

            # sort these now in reverse order, so we favor picking latest links first. This ensures that throughput is grabbed in order
            for lnks_o in lnks_o_by_sat_indx:
                lnks_o.sort(key= lambda l: l.center,reverse=True)

            # for each iteration of this loop, we create a new dr
            while remaining_obs_dv > self.dv_epsilon:
                delta_obs_dv = remaining_obs_dv

                # use center times for now. Later when we validate the data routes we'll check start/end time overlap validity
                curr_time = obs.center
                curr_sat_indx = obs.sat_indx

                # stores windows for dr, in order
                route_winds = [obs]
                route_window_start_sats = {obs: curr_sat_indx}

                iterations = 0
                # keep going till we find a dlnk to end the route
                while True:
                    found_lnk = False
                    for lnk in lnks_o_by_sat_indx[curr_sat_indx]:
                        # find a link that is after current time and has remaining dv routable
                        if lnk.center > curr_time and remaining_dv_o_by_lnk[lnk] > 0:
                            found_lnk = True
                            break

                    # if we have remaining dv to route and we didn't find a link, there's a problem...
                    assert(found_lnk)

                    # we've found our link, move DV along it
                    delta_obs_dv = min(delta_obs_dv,remaining_dv_o_by_lnk[lnk])
                    route_winds.append(lnk)
                    route_window_start_sats[lnk] = curr_sat_indx

                    # print('route_winds')
                    # print(route_winds)

                    if type(lnk) == XlnkWindow:
                        curr_time = lnk.center
                        curr_sat_indx = lnk.rx_sat

                    # if we're at a dlnk, that's the end of the line for this route. Pinch it off as a dr
                    elif type(lnk) == DlnkWindow:
                        remaining_obs_dv -= delta_obs_dv
                        found_routes.append ( 
                            DataRoute(
                                dr_uid, 
                                route =route_winds, 
                                window_start_sats=route_window_start_sats,
                                dv=delta_obs_dv,
                                dv_epsilon =self.dv_epsilon,
                            )
                        )
                        dr_uid += 1

                        # mark that dv as occupied
                        for wind in route_winds:
                            if not type(wind) == ObsWindow:
                                remaining_dv_o_by_lnk[wind] -= delta_obs_dv
                        
                        print('found_routes[-1]')
                        print(found_routes[-1])

                        # from circinus_tools import debug_tools
                        # debug_tools.debug_breakpt()

                        # stop while True loop
                        break

                    iterations += 1
                    if iterations > 100:
                        raise RuntimeWarning('Iterations > 100, likely caught in infinite loop')

            # sanity check that all assigned dv for all routes has been used
            for lnk,remaining_dv in remaining_dv_o_by_lnk.items():
                assert(remaining_dv < self.dv_epsilon)


        for dr in found_routes:
            dr.validate()
            print(str(dr))

        return found_routes, dr_uid



    def extract_utilized_routes( self, copy_routes = True, verbose = False):

        data_routes,dr_uid = self.fabricate_routes()

        from circinus_tools import debug_tools
        debug_tools.debug_breakpt()

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

