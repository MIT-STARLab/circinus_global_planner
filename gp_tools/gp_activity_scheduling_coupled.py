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
from circinus_tools  import time_tools as tt
from circinus_tools  import  constants as const
from circinus_tools.scheduling.custom_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow
from circinus_tools.scheduling.schedule_objects import Dancecard
from circinus_tools.scheduling.routing_objects import DataRoute,DataMultiRoute

def print_verbose(string,verbose=False):
    if verbose:
        print(string)

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

        # number of observations found after filtering
        self.num_obs_filt = const.UNASSIGNED

    def get_stats(self,verbose=True):

        valid = True
        if self.num_obs_filt == 0:
            if verbose:
                print('No observations found! No stats to extract')
            valid = False

        stats = {}
        stats['num_acts'] = len ( self.all_acts_windids) if valid else 0
        stats['num_obs'] = len ( self.all_obs_windids) if valid else 0
        stats['num_dlnks'] = len ( self.all_dlnk_windids) if valid else 0
        stats['num_xlnks'] = len ( self.all_xlnk_windids) if valid else 0
        stats['num_variables'] = self.model.nvariables () if valid else 0
        stats['num_nobjectives'] = self.model.nobjectives () if valid else 0
        stats['num_nconstraints'] = self.model.nconstraints () if valid else 0

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
                if  wind.start >= self.planning_start_dt  and  wind.end  <= self.planning_end_obs_dt:
                    obs_winds_filtered[sat_indx]. append ( wind)

            for wind in dlnk_winds_flat[sat_indx]:
                if  wind.start >= self.planning_start_dt  and  wind.end  <= self.planning_end_dlnk_dt:
                    dlink_winds_flat_filtered[sat_indx]. append ( wind)

            for wind in xlnk_winds_flat[sat_indx]:
                if  wind.start >= self.planning_start_dt  and  wind.end  <= self.planning_end_xlnk_dt:
                    xlink_winds_flat_filtered[sat_indx]. append ( wind)

        return obs_winds_filtered, dlink_winds_flat_filtered, xlink_winds_flat_filtered

    def get_act_model_objs(self,act):

        model_objs_act = {
            'act_object': act,
            'var_dv_utilization': self.model.var_act_dv_utilization[act.window_ID],
            'par_dv_capacity': self.model.par_act_capacity[act.window_ID],
            'var_act_indic': self.model.var_act_indic[act.window_ID],
        }

        return model_objs_act

    def get_dlnk_latency_score_factors(self,all_act_objs_by_windid,obs_windids,dlnk_windids_by_obs_windid,obs_winds_filt,dlnk_winds_filt):

        latency_sf_by_dlnk_obs_windids = {}

        # loop through all obs - dlnk combinations and determine latency score factors for them
        for o in obs_windids:
            latency_by_dlnk_windid = {}

            obs = all_act_objs_by_windid[o]
            for d in dlnk_windids_by_obs_windid[o]:
                dlnk = all_act_objs_by_windid[d]
                latency_by_dlnk_windid[d] = DataRoute.calc_latency(
                        obs,
                        dlnk,
                        units='minutes',
                        obs_option = self.latency_params['obs'], 
                        dlnk_option = self.latency_params['dlnk']
                    )

            #  the shortest latency downlink for this observation has a score factor of 1.0, and the score factors for the other downlinks decrease as the inverse of increasing latency
            min_lat = min(latency_by_dlnk_windid.values())
            min_lat = max(min_lat, self.min_latency_for_sf_1_mins)
            for d,lat in latency_by_dlnk_windid.items():
                lat = max(lat, self.min_latency_for_sf_1_mins)
                latency_sf_by_dlnk_obs_windids[(d,o)] = min_lat/lat
 
        return latency_sf_by_dlnk_obs_windids

    def get_activity_structs(self,obs_winds_filt,dlnk_winds_filt,xlnk_winds_filt):
        """Get the structures representing activity subscripts for variables, parameters that are used in creating the pyomo model"""

        #  all activities are uniquely indexed. these structures keep track of those, and the mapping to activity objects
        all_acts_windids = set()
        all_obs_windids = set()
        all_dlnk_windids = set()
        all_xlnk_windids = set()
        all_act_objs_by_windid = {}
        all_act_windids_by_obj = {}
        capacity_by_act_windid = {}

        # all acts for a given sat
        sats_acts_set = [set() for sat_indx in range (self.num_sats)]
        sats_dlnks_set = [set() for sat_indx in range (self.num_sats)]

        # this is T^s - defines a unique time system for every satellite based on where the activities actually take place, as opposed to a finely granularized, monolithic list of times that could impose MILP constraints even when there are no activities happening during most of the times.
        all_link_times_by_sat_indx = {sat_indx: [] for sat_indx in range(self.num_sats)}  
        # this is used for producing T^s
        all_link_times_link_objs_by_sat_indx = {sat_indx: [] for sat_indx in range(self.num_sats)}  
        # this is a list of all the (s,o,t) in v_(s,o,t) in the formulation (where t is an index, starting at 0 for each sat)
        sat_obs_time_subscripts = []
        # same as above, but just (s,t)
        sat_time_subscripts = set()
        # this is the actual time corresponding to the subscript stored in sat_obs_time_subscripts
        time_by_sat_time_subscript = {}
        # keeps track of the v_(s,o,t) subscripts corresponding to dv availability for a given activity (v_(x,o)^-, v_(d,o))
        sat_obs_time_subscripts_by_link_obj = {}
        # opposite direction of the above - look up link objects by v_(s,o,t) subscript
        link_objs_by_sat_obs_time_subscripts = {}
        # all the obs found at (s,t)
        obs_windids_by_sat_time_subscript = {}

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
                all_act_windids_by_obj[dlnk] = dlnk_windid
                all_acts_windids.add(dlnk_windid)
                all_dlnk_windids.add(dlnk_windid)
                capacity_by_act_windid[dlnk_windid] = dlnk.data_vol
                obs_windids_by_dlnk_windid.setdefault(dlnk_windid, set())
                all_link_times_link_objs_by_sat_indx[sat_indx].append((dlnk.center,[dlnk]))
                sat_obs_time_subscripts_by_link_obj.setdefault(dlnk, [])
                sats_acts_set[sat_indx].add(dlnk)
                sats_dlnks_set[sat_indx].add(dlnk)

                for obs_sat in obs_winds_filt:
                    for obs in obs_sat:
                        # we'll do this obs storing again for every subsequent downlink, of course...it's fine, we'll survive
                        obs_windid = obs.window_ID
                        all_act_objs_by_windid[obs_windid] = obs
                        all_act_windids_by_obj[obs] = obs_windid
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
                all_act_windids_by_obj[xlnk] = xlnk_windid
                all_acts_windids.add(xlnk_windid)
                all_xlnk_windids.add(xlnk_windid)
                capacity_by_act_windid[xlnk_windid] = xlnk.data_vol
                # need set() here because xlnk windows are seen by sat_indcs on both sides of the link
                obs_windids_by_xlnk_windid.setdefault(xlnk_windid, set())
                # note that we can't use xlnk.sat_indx below below every xlnk wind has a sat_indx and an xsat_indx (sat_indx < xsat_indx) and it's incorrect to use sat_indx straight
                all_link_times_link_objs_by_sat_indx[sat_indx].append((xlnk.center,[xlnk]))
                sats_acts_set[sat_indx].add(xlnk)
                sat_obs_time_subscripts_by_link_obj.setdefault(xlnk, [])

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

        # deal with the fact that some activities have pretty much the same center time, and we need to merge those into a single timepoint in T_s rather than multiple. This is necessary because we need to be able to sum together the effects of all links that are overlapping in time, so that we can account for their combined effect in constraints. E.g., the sum of data volume sent out from a satellite over two outgoing, concurrent crosslinks <= available dv before their mutual center time. Note that this is not strictly necessary for when we constrain each satellite to only allow a single activity at a time, but it's good to build this in now so that the underlying model is safer and less able to bug out.
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

        # generate list of subscripts for v_(s,o,t), as well as the correspondence between activity windows and their entries in v_(s,o,t)
        for obs_sat in obs_winds_filt:
            for obs in obs_sat:
                for sat_indx in range(self.num_sats):
                    for time_indx,(time,link_objs) in enumerate(all_link_times_link_objs_by_sat_indx_merged[sat_indx]):
                        sat_time_subscr = (sat_indx,time_indx)
                        time_by_sat_time_subscript[sat_time_subscr] = time
                        obs_windids_by_sat_time_subscript.setdefault(sat_time_subscr,[])

                        if time > obs.center:
                            sat_obs_time_subscr = (sat_indx,obs.window_ID,time_indx)
                            sat_obs_time_subscripts.append(sat_obs_time_subscr)
                            sat_time_subscripts.add(sat_time_subscr)
                            obs_windids_by_sat_time_subscript[sat_time_subscr].append(obs.window_ID)
                            link_objs_by_sat_obs_time_subscripts[sat_obs_time_subscr] = link_objs
                            
                            for link_obj in link_objs:
                                # if sat_indx ==4 and obs.window_ID == 41 and link_obj.window_ID == 5303:
                                #     from circinus_tools import debug_tools
                                #     debug_tools.debug_breakpt()

                                sat_obs_time_subscripts_by_link_obj[link_obj].append(sat_obs_time_subscr)


        sats_acts = []
        sats_dlnks = []
        for sat_indx in range(self.num_sats):

            sat_act_list = list(sats_acts_set[sat_indx])
            sat_act_list.sort(key= lambda a: a.center)
            sats_acts.append(sat_act_list)

            sat_dlnk_list = list(sats_dlnks_set[sat_indx])
            sat_dlnk_list.sort(key= lambda a: a.center)
            sats_dlnks.append(sat_dlnk_list)

        return sats_acts,sats_dlnks,all_acts_windids,all_obs_windids,all_dlnk_windids,all_xlnk_windids,all_act_objs_by_windid,all_act_windids_by_obj,capacity_by_act_windid,dlnk_obs_windids,xlnk_obs_windids,obs_windids_by_lnk_windid,obs_windids_by_dlnk_windid,obs_windids_by_xlnk_windid,dlnk_windids_by_obs_windid,xlnk_windids_by_obs_windid,all_link_times_by_sat_indx,sat_obs_time_subscripts,time_by_sat_time_subscript,sat_obs_time_subscripts_by_link_obj,link_objs_by_sat_obs_time_subscripts,sat_time_subscripts,obs_windids_by_sat_time_subscript


    def make_model ( self,obs_winds,dlnk_winds_flat,xlnk_winds_flat, ecl_winds, verbose = False):
        # note this model only works for non-symmetric crosslink windows!

        model = pe.ConcreteModel()
        self.model = model

        # filter the routes to make sure that  none of their activities fall outside the scheduling window
        obs_winds_filt,dlnk_winds_filt,xlnk_winds_filt = self.filter_windows(obs_winds,dlnk_winds_flat,xlnk_winds_flat,self.num_sats)

        self.ecl_winds = ecl_winds

        self.num_obs_filt = sum(len(obs_winds_sat) for obs_winds_sat in obs_winds_filt)
        if self.num_obs_filt == 0:
            if verbose:
                print('No observations found! Quitting coupled AS early')
            self.model_constructed = False
            return

        print_verbose('num obs winds filt: %d'%(sum(len(winds) for winds in obs_winds_filt)),verbose)
        print_verbose('num dlnk winds filt: %d'%(sum(len(winds) for winds in dlnk_winds_filt)),verbose)
        print_verbose('num xlnk winds filt: %d'%(sum(len(winds) for winds in xlnk_winds_filt)),verbose)

        ##############################
        #  Make indices/ subscripts
        ##############################

        try:
            (sats_acts,
                sats_dlnks,
                all_acts_windids,
                all_obs_windids,
                all_dlnk_windids,
                all_xlnk_windids,
                all_act_objs_by_windid,
                all_act_windids_by_obj,
                capacity_by_act_windid,
                dlnk_obs_windids,
                xlnk_obs_windids,
                obs_windids_by_lnk_windid,
                obs_windids_by_dlnk_windid,
                obs_windids_by_xlnk_windid,
                dlnk_windids_by_obs_windid,
                xlnk_windids_by_obs_windid,
                all_link_times_by_sat_indx,
                sat_obs_time_subscripts,
                time_by_sat_time_subscript,
                sat_obs_time_subscripts_by_link_obj,
                link_objs_by_sat_obs_time_subscripts,
                sat_time_subscripts,
                obs_windids_by_sat_time_subscript
                 ) = self.get_activity_structs(
                        obs_winds_filt,
                        dlnk_winds_filt,
                        xlnk_winds_filt
                    )

            # verify that all acts found are within the planning window, otherwise we may end up with strange results
            for sat_acts in sats_acts:
                for act in sat_acts:
                    if act.start < self.planning_start_dt or act.end > self.planning_end_dt:
                        raise RuntimeWarning('Activity is out of planning window range (start %s, end %s): %s'%(self.planning_start_dt,self.planning_end_dt,act))

            self.all_acts_windids = all_acts_windids
            self.all_obs_windids = all_obs_windids
            self.all_dlnk_windids = all_acts_windids
            self.all_xlnk_windids = all_acts_windids
            self.all_act_objs_by_windid = all_act_objs_by_windid
            self.time_by_sat_time_subscript = time_by_sat_time_subscript
            self.all_link_times_by_sat_indx = all_link_times_by_sat_indx
            self.obs_windids_by_sat_time_subscript = obs_windids_by_sat_time_subscript
            # combine xlnk and dlnks windids by obs windid together
            self.lnk_windids_by_obs_windid = copy(xlnk_windids_by_obs_windid)
            for obs_windid,dlnk_windids in dlnk_windids_by_obs_windid.items():
                self.lnk_windids_by_obs_windid.setdefault(obs_windid, set())
                self.lnk_windids_by_obs_windid[obs_windid] = self.lnk_windids_by_obs_windid[obs_windid].union(dlnk_windids)


            latency_sf_by_dlnk_obs_windids =  self.get_dlnk_latency_score_factors(
                all_act_objs_by_windid,
                all_obs_windids,
                dlnk_windids_by_obs_windid,
                obs_winds_filt,
                dlnk_winds_filt
            )

            # construct a set of dance cards for every satellite, 
            # each of which keeps track of all of the activities of satellite 
            # can possibly execute at any given time slice delta T. 
            # this is for constructing energy storage constraints
            # using resource_delta_t_s because this dancecard is solely for use in constructing resource constraints
            # note that these dancecards will baloon in size pretty quickly as the planning_end_dt increases. However most of the complexity they introduce is before planning_end_obs,xlnk_dt because that's the horizon where obs,xlnk actitivities are included. After that there should only be sparse downlinks
            es_act_dancecards = [Dancecard(self.planning_start_dt,self.planning_end_dt,self.resource_delta_t_s,item_init=None,mode='timestep') for sat_indx in range (self.num_sats)]
            
            for sat_indx in range (self.num_sats): 
                es_act_dancecards[sat_indx].add_winds_to_dancecard(sats_acts[sat_indx])
                es_act_dancecards[sat_indx].add_winds_to_dancecard(ecl_winds[sat_indx])


        except IndexError:
            raise RuntimeWarning('sat_indx out of range. Are you sure all of your input files are consistent? (including pickles)')        
        #  subscript for each activity a
        model.act_windids = pe.Set(initialize= all_acts_windids)
        #  subscript for each obs o
        model.obs_windids = pe.Set(initialize= all_obs_windids)
        #  subscript for each link
        model.lnk_windids = pe.Set(initialize= all_dlnk_windids.union(all_xlnk_windids))
        # every act-obs within that act combination
        model.lnk_obs_subscripts = pe.Set(initialize= dlnk_obs_windids.union(xlnk_obs_windids),dimen=2)
        # every possible dlnk-obs combination
        model.dlnk_obs_subscripts = pe.Set(initialize= dlnk_obs_windids,dimen=2)
        #  subscript for each satellite
        model.sat_indcs = pe.Set(initialize= range(self.num_sats))
        #  subscript for each s,o,t in v_(s,o,t)
        # t in model.sat_obs_time_subscripts is index into all_link_times_by_sat_indx
        model.sat_obs_time_subscripts = pe.Set(initialize= sat_obs_time_subscripts,dimen=3)


        # timepoints is the indices, which starts at 0 
        #  NOTE: we assume the same time system for every satellite
        self.es_time_getter_dc = es_act_dancecards[0]
        es_num_timepoints = es_act_dancecards[0].num_timepoints
        model.es_timepoint_indcs = pe.Set(initialize=  self.es_time_getter_dc.get_tp_indcs())

        if self.solver_name == 'gurobi' or self.solver_name == 'cplex':
            int_feas_tol = self.solver_params['integer_feasibility_tolerance']
        else:
            raise NotImplementedError

        for d_o,lat_sf in latency_sf_by_dlnk_obs_windids.items():        
            if lat_sf > int_feas_tol*self.big_M_lat:
                raise RuntimeWarning('big_M_lat (%f) is not large enough for latency score factor %f and integer feasibility tolerance %f (dlnk windid %d obs windid %d)'%(self.big_M_lat,lat_sf,int_feas_tol,d_o[0],d_o[1]))


        ##############################
        #  Make parameters
        ##############################

        model.par_min_obs_dv_dlnk_req = pe.Param (initialize=self.min_obs_dv_dlnk_req)
        model.par_total_obs_dv = sum(capacity_by_act_windid[o] for o in all_obs_windids)
        model.par_act_capacity = pe.Param(model.act_windids,initialize =capacity_by_act_windid)

        model.par_latency_sf_dlnk_obs = pe.Param(model.dlnk_obs_subscripts,initialize=latency_sf_by_dlnk_obs_windids)

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

        # v_a, variables [1]
        model.var_act_dv_utilization  = pe.Var (model.act_windids, within = pe.NonNegativeReals)
        # v_(a,o), variables [2]
        model.var_lnk_obs_dv_utilization  = pe.Var (model.lnk_obs_subscripts,  within = pe.NonNegativeReals)
        # v_(s,o,t), variables [3]
        model.var_sat_obs_time_dv  = pe.Var (model.sat_obs_time_subscripts,  within = pe.NonNegativeReals)
        # I_a, variables [6]
        model.var_act_indic  = pe.Var (model.act_windids, within = pe.Binary)
        # satellite energy storage, e_(s,t), variables [7]
        model.var_sats_estore  = pe.Var (model.sat_indcs,  model.es_timepoint_indcs,  within = pe.NonNegativeReals)
        # I_(d,o), variables [5]
        model.var_dlnk_obs_indic  = pe.Var (model.dlnk_obs_subscripts, within = pe.Binary)
        # f_o, variables [4]
        model.var_latency_sf_obs = pe.Var (model.obs_windids,  bounds = (0,1.0))
        
        if self.allow_act_timing_constr_violations:
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

        print_verbose('creating constraints',verbose)

        # TODO: renumber  these with the final numbering

        # note that the observations show up within model.acts as well, so we also constraint route scheduled DV by the real available DV from each observation
        def c1_rule( model,a):
            return ( model.var_act_dv_utilization[a] <= model.par_act_capacity[a])
        model.c1 =pe.Constraint ( model.act_windids,  rule=c1_rule)

        def c17_rule( model,a):
            return model.var_act_indic[a] >=  model.var_act_dv_utilization[a]/model.par_act_capacity[a]
        model.c17 =pe.Constraint ( model.act_windids,  rule=c17_rule) 

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

        print_verbose('make c5 constraints',verbose)

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

                    # note this needs to be an equality in order to use var_sat_obs_time_dv to produce data storage output in extract_resource_usage below
                    model.c5.add(model.var_sat_obs_time_dv[s,o,t_indx] == dv_sum)


        # c6
        model.c6  = pe.ConstraintList()

        for sat_obs_time_subscript,link_objs in link_objs_by_sat_obs_time_subscripts.items():
            # sum of all the links with outgoing dv at s,o,t
            outgoing_dv_sum = 0

            s = sat_obs_time_subscript[0]
            o = sat_obs_time_subscript[1]

            for lnk in link_objs:
                # verify lnk is actually a link act
                if not (type(lnk) == DlnkWindow or type(lnk) == XlnkWindow):
                    continue

                # if receiving, we're not constrained by any dv availability
                # need to consider this here because the link_objs_by_sat_obs_time_subscripts mapping holds both tx and rx directions for every link
                if type(lnk) == XlnkWindow: 
                    # this code won't work if we're assuming symmetric xlnk windows
                    assert(not lnk.symmetric)
                    if lnk.is_rx(s):
                        continue # corresponds to for lnk in link_objs loop

                outgoing_dv_sum += model.var_lnk_obs_dv_utilization[lnk.window_ID,o]

            model.c6.add( outgoing_dv_sum <= model.var_sat_obs_time_dv[sat_obs_time_subscript])


        def dlnk_lat_time_getter(dlnk):
            return getattr(dlnk,self.latency_params['dlnk'])

        #  obs downlink indicator constraints [8]
        model.c8  = pe.ConstraintList()
        for o in model.obs_windids:
            dlnk_windids_sorted = sorted(dlnk_windids_by_obs_windid[o],key = lambda d: dlnk_lat_time_getter(all_act_objs_by_windid[d]))

            for dlnk_indx,d in enumerate(dlnk_windids_sorted):
                dlnks_cum_dv_up_to = sum(model.var_lnk_obs_dv_utilization[dlnk_windid_up_to,o] for dlnk_windid_up_to in dlnk_windids_sorted[:(dlnk_indx+1)])

                model.c8.add( dlnks_cum_dv_up_to >= model.par_min_obs_dv_dlnk_req * model.var_dlnk_obs_indic[d,o])




        #  observation latency score factor constraints [9]
        model.c9  = pe.ConstraintList()
        for o in model.obs_windids:

            # sort the latency score factors for all the potential downlinks for this observation in increasing order -  important for constraint construction
            dlnk_windids_sorted = sorted(dlnk_windids_by_obs_windid[o],key = lambda d: model.par_latency_sf_dlnk_obs[d,o])

            #  initial constraint -  score factor for this observation will be equal to zero if no dmrs for this obs were chosen
            model.c9.add( model.var_latency_sf_obs[o] <= 0 + self.big_M_lat * sum(model.var_dlnk_obs_indic[d,o] for d in dlnk_windids_sorted))
            
            for curr_d_indx,curr_d in enumerate(dlnk_windids_sorted):
                #  add constraint that score factor for observation is less than or equal to the score factor for this dlnk, plus any big M terms for any dlnks with larger score factors.
                #  what this does is effectively disable the constraint for the score factor for this dlnk if any higher score factor dlnks was chosen to fulfill minimum dv downlink requirement 
                model.c9.add( model.var_latency_sf_obs[o] <= 
                    model.par_latency_sf_dlnk_obs[curr_d,o] + 
                    self.big_M_lat * sum(model.var_dlnk_obs_indic[d,o] for d in dlnk_windids_sorted[curr_d_indx+1:]) )

        print_verbose('make act overlap constraints',verbose)

        #  intra-satellite activity overlap constraints [10],[11],[12]
        #  well, 12 is activity minimum time duration
        model.c10_11  = pe.ConstraintList()
        model.c12  = pe.ConstraintList()
        # pass the model objects getter function so it can be called in place
        self.gen_intra_sat_act_overlap_constraints(model.c10_11,model.c12,sats_acts,self.get_act_model_objs,constraint_violation_model_objs)

        # inter-satellite downlink overlap constraints [9],[10]
        model.c14_15  = pe.ConstraintList()
        # pass the model objects getter function so it can be called in place
        self.gen_inter_sat_act_overlap_constraints(model.c14_15,sats_dlnks,self.get_act_model_objs,constraint_violation_model_objs)

        #  energy constraints [13]
        # todo: maybe this ought to be moved to the super class, but i don't anticipate this code changing much any time soon, so i'll punt that.
        model.c13  = pe.ConstraintList()
        for sat_indx in range (self.num_sats): 

            # tp_indx serves as an index into the satellite activity dance cards
            for tp_indx in model.es_timepoint_indcs:
                #  constraining first time step to initial energy storage
                #  continue for loop afterwards because no geq/leq constraints needed for this index
                if tp_indx == 0:
                    model.c13.add( model.var_sats_estore[sat_indx,tp_indx] ==  model.par_sats_estore_initial[sat_indx])
                    continue 

                #  minimum and maximum storage constraints
                model.c13.add( model.var_sats_estore[sat_indx,tp_indx] >= model.par_sats_estore_min[sat_indx])
                model.c13.add( model.var_sats_estore[sat_indx,tp_indx] <= model.par_sats_estore_max[sat_indx])

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
                                    * model.var_act_dv_utilization[act_windid]/model.par_act_capacity[act_windid]
                                    * model.par_resource_delta_t
                                )

                            #  if the satellite is not in sunlight then we can't charge
                            elif type(act) == EclipseWindow:
                                charging = False

                    # add in charging energy contribution ( if possible)
                    charging_delta_e = model.par_sats_edot_by_mode[sat_indx]['orbit_insunlight_average_charging']*model.par_resource_delta_t if charging else 0

                    #  base-level satellite energy usage (not including additional activities)
                    base_delta_e = model.par_sats_edot_by_mode[sat_indx]['base']*model.par_resource_delta_t

                    # maximum bound of energy at current time step based on last time step
                    model.c13.add( model.var_sats_estore[sat_indx,tp_indx] <= 
                        model.var_sats_estore[sat_indx,tp_indx-1]
                        + activity_delta_e
                        + charging_delta_e
                        + base_delta_e
                    )

                    # maximum bound of energy at current time step based on last time step
                    model.c13.add( model.var_sats_estore[sat_indx,tp_indx] >= 
                        model.var_sats_estore[sat_indx,tp_indx-1]
                        + activity_delta_e
                        + base_delta_e
                    )


        #  data storage constraints [16]
        model.c16  = pe.ConstraintList()
        if self.enforce_data_storage_constr:
            for sat_time_subscr in sat_time_subscripts:
                dv_sum = 0
                for obs_windid in obs_windids_by_sat_time_subscript[sat_time_subscr]:
                    sat_indx = sat_time_subscr[0]
                    time_indx = sat_time_subscr[1]
                    sat_obs_time_subscr = (sat_indx,obs_windid,time_indx)
                    dv_sum += model.var_sat_obs_time_dv[sat_obs_time_subscr]
                
                model.c16.add( dv_sum <= model.par_sats_dstore_max[sat_indx])



        # note: c17 is way above 



        #  note: use model.c8[indx].expr.to_string()  to print out the constraint in a human readable form
        #                ^ USES BASE 1 INDEXING!!! WTF??
                

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        ##############################
        #  Make objective
        ##############################

        print_verbose('make obj',verbose)

        #  determine which time points to use as "spot checks" on resource margin. These are the points that will be used in the objective function for maximizing resource margin
        timepoint_spacing = ceil(es_num_timepoints/self.resource_margin_obj_num_timepoints)
        # need to turn the generator into a list for slicing
        #  note: have to get the generator again
        decimated_tp_indcs = list(self.es_time_getter_dc.get_tp_indcs())[::timepoint_spacing]
        rsrc_norm_f = len(decimated_tp_indcs) * len(model.sat_indcs)

        def obj_rule(model):
            total_dv_term = self.obj_weights['obs_dv'] * 1/model.par_total_obs_dv * sum(model.var_act_dv_utilization[o] for o in model.obs_windids) 

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

            return total_dv_term + latency_term + energy_margin_term - inter_sat_act_constr_violations_term - intra_sat_act_constr_violations_term
            # return total_dv_term + latency_term + energy_margin_term
            
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

            if str (v) =='var_dlnk_obs_indic': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" %s %d   %s %s"%(index, val, self.all_act_objs_by_windid[index[0]],self.all_act_objs_by_windid[index[1]]))


            print_times = False
            if print_times:
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
                        #     print (" %s %0.3f   %s"%(index, val, self.time_by_sat_time_subscript[index]))

                    for sat_indx in range(self.num_sats):
                        # sort by time index
                        sorted_var_sat_obs_time_dv[sat_indx].sort(key= lambda subscr: subscr[2])

                        for index in sorted_var_sat_obs_time_dv[sat_indx]:
                            val  = varobject[index].value
                            obs_windid = index[1]
                            time_indx = index[2]
                            # if sat_indx == 4 and  obs_windid == 41:
                            if sat_indx == 4 and obs_windid == 41:
                                print (" %s %0.3f   %s"%(index, val, self.time_by_sat_time_subscript[(sat_indx,time_indx)]))

    def fabricate_routes(self):

        found_routes = []
        dr_uid = 0
        for o in self.model.obs_windids:

            # if not enough dv from this obs, then don't construct routes for it
            if not (pe.value(self.model.var_act_dv_utilization[o]) > self.min_obs_dv_dlnk_req):
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

            # sort these now in order, so we favor picking earliest links first. This ensures that we choose to construct data routes through the hardest-to-reach (earliest) links first
            for lnks_o in lnks_o_by_sat_indx:
                lnks_o.sort(key= lambda l: l.center)

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
                        if lnk.center > curr_time and remaining_dv_o_by_lnk[lnk] > self.dv_epsilon:
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
                                self.gp_agent_ID,
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
                        
                        # stop while True loop
                        break

                    iterations += 1
                    if iterations > 100:
                        raise RuntimeWarning('Iterations > 100, likely caught in infinite loop')

            # sanity check that all assigned dv for all routes has been used
            for lnk,remaining_dv in remaining_dv_o_by_lnk.items():
                assert(remaining_dv < self.dv_epsilon)


        # only validate based on center time for now
        for dr in found_routes:
            dr.validate(time_option='center')
            # print(str(dr))

        return found_routes, dr_uid



    def extract_utilized_routes( self, copy_routes = True, verbose = False):

        if verbose:
            print ('utilized routes:')

        scheduled_routes_flat = []

        if self.num_obs_filt == 0:
            if verbose:
                print('No observations found! No routes to extract')
            return scheduled_routes_flat

        # create the basic set of data routes from scheduled windows, dv usage
        data_routes,dr_uid = self.fabricate_routes()

        # commenting out for now because I'm not sure this is actually useful, considering we already have to modify the windows below...
        # def copy_choice(route):
        #     if copy_routes:
        #         # this will copy everything but the underlying activity windows. 
        #         return copy(route)
        #     else:
        #         return route


        # the rest of the code expects DataMultiRoute objects, so we'll create these as shallow wrappers around the data routes we just fabricated
        for dr in data_routes:
            scheduled_route = DataMultiRoute(self.gp_agent_ID,dr_uid,[dr],dv_epsilon=self.dv_epsilon)
            # the dr within is fully utilized (1.0) 
            scheduled_route.set_scheduled_dv_frac(1.0)
            scheduled_routes_flat.append(scheduled_route)
            if verbose:
                print(scheduled_route)
                
            dr_uid+=1

        # the below check code is adapted from gp_activity_scheduling_separate, and might be a bit superfluous...
        # todo: update the code?

        # examine the schedulable data volume for every activity window, checking as we go that the data volume is sufficient for at least the route in which the window is found
        #  note that this code is slightly inefficient because it might duplicate windows across routes. that's fine though, because we're thorough in checking across all routes
        # note: dmr is for DataMultiRoute
        wind_sched_dv_check = {}
        for dmr in scheduled_routes_flat:
            # wind may get set multiple times due to Windows appearing across routes, but that's not really a big deal
            for wind in dmr.get_winds():
                # var_act_dv_utilization is a lower bound on how much data volume was used by the window
                wind_sched_dv_check[wind] = pe.value(self.model.var_act_dv_utilization[wind.window_ID])
                #  initialize this while we're here
                wind.scheduled_data_vol = 0

            # similar to winds, set dr data vols to 0
            for dr in dmr.data_routes:
                # dmrs should not be sharing drs across themselves, so no scheduled dv should be seen yet
                assert(dr.scheduled_dv == const.UNASSIGNED)
                dr.scheduled_dv = dmr.scheduled_dv_by_dr[dr] 


        #  now we want to mark the real scheduled data volume for every window. we only actually need to use as much data volume as the data routes want to push through the window
        #  add data volume for every route passing through every window
        for dmr in scheduled_routes_flat:
            for wind in dmr.get_winds():
                wind.scheduled_data_vol += dmr.scheduled_dv_for_wind(wind)


        # update the window beginning and end times based upon their amount of scheduled data volume
        # keep track of which ones we've updated, because we should only update once
        updated_winds = set()
        for dmr in scheduled_routes_flat:
            for wind in dmr.get_winds():
                #  this check should be at least as big as the scheduled data volume as calculated from all of the route data volumes. (it's not constrained from above, so it could be bigger)
                if wind_sched_dv_check[wind] < wind.scheduled_data_vol - self.dv_epsilon:
                    raise RuntimeWarning('inconsistent activity scheduling results: activity window data volume determined from route dvs (%f) is greater than dv from var_act_indic (%f) [dmr: %s, wind: %s]'%(wind.scheduled_data_vol,wind_sched_dv_check[wind],dmr,wind))

                if not wind in updated_winds:
                    # note that the line below seems like it may break the scheduled times for activities by specifying a minimum activity duration. however, this minimum activity duration is already accounted for in scheduling constraints
                    wind.update_duration_from_scheduled_dv (min_duration_s=self.min_act_duration_s[type(wind)])
                    updated_winds.add(wind)

            # validate the data multi route (and in turn, the scheduled data vols of all the data routes under it)
            if self.allow_act_timing_constr_violations:
                dmr.validate(time_option='center') # this is bad to use in general, allows window overlaps to occur
            else:
                dmr.validate()

        return scheduled_routes_flat

    def extract_resource_usage( self, decimation_factor =1, verbose = False):

        if self.num_obs_filt == 0:
            if verbose:
                print('No observations found! No resource usage data to extract')
            return None,None

        energy_usage = {}

        # note t_vals is relative time to start
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

        # note t_vals is relative time to start
        t_vals = [[] for sat_indx in range ( self.num_sats)]
        d_vals = [[] for sat_indx in range ( self.num_sats)]

        t_val_initial = 0
        # starts out with zero DV todo: fix this if incorrect...
        d_val_initial = 0

        # TODO: this code feels super inefficient somehow.  make it better?
        for sat_indx, sat in enumerate (self.model.sat_indcs):

            t_vals[sat_indx].append(t_val_initial)
            d_vals[sat_indx].append(d_val_initial)
            last_d_val = d_val_initial

            # last_tp_indx = 0
            # for tp_indx in self.model.ds_timepoint_indcs:
            for time_indx,time in enumerate(self.all_link_times_by_sat_indx[sat_indx]):

                # assuming output in minutes for now
                curr_time = (time-self.planning_start_dt).total_seconds()/60

                curr_dv = 0
                for obs_windid in self.obs_windids_by_sat_time_subscript[(sat_indx,time_indx)]:
                    curr_dv += pe.value(self.model.var_sat_obs_time_dv[(sat_indx,obs_windid,time_indx)])


                # note that these times and dvs are based on activity center times (we're repurposing the variables used in constraint construction, for output). This is a bit of an artifact, but is okay because we DO meet DV constraints. This will show up on the output plot as an inprecise match between the times of executed activities and the changes in stored DV they cause. This is because we constrain the DV used for a link at timepoint t by the DV contributions from all acts up to t-1. That means that the DV contribution of a crosslink, say, that was in the middle of executing at t-1 will only be seen at t.
                # This is similar to the situation in the separate AS model
                # todo: if a more precise accounting of the dv on each satellite at each time is needed, then one should go through all windows and do bookkeeping on when exactly DV is moved from sat to sat.

                # want the output to be represented as a step change, so add the last dv val, at the current time, as well
                t_vals[sat_indx].append(curr_time)
                d_vals[sat_indx].append(last_d_val)
                t_vals[sat_indx].append(curr_time)
                d_vals[sat_indx].append(curr_dv)

                last_d_val = curr_dv

            # also add last time
            curr_time = (self.planning_end_dt-self.planning_start_dt).total_seconds()/60
            t_vals[sat_indx].append(curr_time)
            d_vals[sat_indx].append(last_d_val)


        data_usage['time_mins'] = t_vals
        data_usage['d_sats'] = d_vals

        return  energy_usage, data_usage
