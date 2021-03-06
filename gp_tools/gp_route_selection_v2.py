# Algorithm for creating a set of data routes for a given observation, later fed into the activity scheduling stage
# 
# This version uses a dynamic programming algorithm to optimize selected paths
# 
# @author Kit Kennedy
#
#  note that a path is the same as a route. 

from  datetime import timedelta
from copy import copy,deepcopy
from collections import namedtuple
from numpy import argmin
from math import floor

from scipy.optimize import linprog

# for Birdseye use
# from cmdline:
# $ birdseye
# http://localhost:7777/
# in source:
# from birdseye import eye
# @eye
# def yo_mama():
# ...

from circinus_tools  import time_tools as tt
from circinus_tools.scheduling.custom_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow
from circinus_tools.scheduling.schedule_objects import Dancecard
from circinus_tools.scheduling.routing_objects import DataRoute, DataMultiRoute
from circinus_tools.scheduling.base_window import get_pairwise_overlap_max_dv
from circinus_tools.activity_bespoke_handling import ActivityTimingHelper

from circinus_tools import debug_tools


class DeconflictedRoute():
    """docstring for DeconflictedRoute"""
    def __init__(self,dr,available_dv):
        self.dr = dr
        self.available_dv = available_dv

    def __repr__( self):
        return 'DeconflictedRoute(%s,%s)'%( self.dr, self.available_dv)

class RouteRecord():
    """docstring for RouteRecord"""
    def __init__(self,dv,release_time,routes=[]):
        """create a new route record
        
        used in dynamic programming algorithm to keep track of data volume and data routes leading to a given point in time on a satellite
        :param dv:  data volume for set of routes
        :type dv:  float
        :param routes:  list of DataRoute objects that currently end at a given time on a satellite , defaults to []
        :type routes: list, optional
        """
        self.dv = dv
        # self._routes is sorted by route end time (end time of last act in route), in increasing order
        self._routes = routes
        self.release_time = release_time

        self.output_date_str_format = 'short'


    def __repr__( self):
        return 'RouteRecord(%s,%s,%s)'%( self.dv, tt.date_string(self.release_time,self.output_date_str_format),self._routes)

    def __copy__(self):
        newone = type(self)( self.dv,copy(self.release_time))
        #  make a shallow copy of every data route object - this avoids copying the underlying windows image data route
        newone._routes = [copy(rt) for rt in self._routes]
        return newone

    def __iter__(self):
        return (dr for dr in self._routes)

    def __len__(self):
        return len(self._routes)

    def add_to_routes(self,new_dr):
        """Add to the routes stored in this route record, sorting if need be"""

        latest_dr_end = None
        if len(self._routes) > 0:
            latest_dr_end = self._routes[-1].get_end().original_end

        end_new_dr = new_dr.get_end().original_end

        self._routes.append(new_dr)

        # sort by end time of the routes
        # ( only if we have to, because it's expensive)
        if latest_dr_end is not None and end_new_dr < latest_dr_end:
            self._routes.sort(key=lambda dr:dr.get_end().original_end)

    def get_act_dv_availability(self,act,sat_indx,act_timing_helper):
        act_start = act.original_start
        act_center = act.center

        # look through routes in reverse order of end time
        # self._routes needs to be sorted so that returning False without looking through whole list is still accurate.
        # for dr in reversed(self._routes):

        dv_avail = act.data_vol

        reduce_dv_for_overlap = True

        for dr in self._routes:
            # todo: assess if this code takes way too long because it checks every dr in self._routes

            dr_end_act = dr.get_end()

            transition_time_req = act_timing_helper.get_transition_time_req(dr_end_act,act,sat_indx,sat_indx)

            # don't calculate overlapped dv if there's sufficient transition time
            if (act_start - dr_end_act.original_end).total_seconds() >= transition_time_req:
                continue
            elif (act.center - dr_end_act.center).total_seconds() < transition_time_req:
                dv_avail = 0
                continue

            # if we want to actually reduce activity data volume (as marked in data route that will be created) to reflect overlap. Otherwise just return the activity dv.
            # it's a good policy to reduce the dv avail to reflect the transition time (to not make routes seem more promising than they actually are), but can take a while to calculate
            if reduce_dv_for_overlap:
                dv_avail_with_dr = get_pairwise_overlap_max_dv(dr_end_act, act, transition_time_req)

                dv_avail = min(dv_avail,dv_avail_with_dr)


        return dv_avail

    def get_dv_deconflicted_routes(self,rr_other,routable_obs_dv_multiplier = 3, min_dv=0,verbose = False):
        """ determines which routes from other route record are able to be extended to the calling route record, from a data volume availability standpoint.
        
        determines from the set of routes contained in self which routes in rr_other are able to be extended to the satellite and time point for the self route record. determines this by figuring out what data volume is still available in windows to run routes through, and optimally selecting which routes from rr_other are allocated that data volume. optimal here means the allocation that maximizes the data volume available in the de-conflicted routes returned. Note that this function DOES NOT currently look at transition time conflicts - temporal overlap is not checked, only if the routes share the same window object or not.

        :param rr_other:  the other route record, containing a set of data routes possible to extend to self
        :type rr_other: RouteRecord
        :param min_dv:  minimum data volume allowable for a route to be selectable, defaults to 0
        :type min_dv: number, optional
        :param verbose: output verbosity, defaults to False
        :type verbose: bool, optional
        :returns:  a set of de-conflicted routes with data volume allocations that, if all selected together, would violate no constraints on window data volume capacities across all routes in self
        :rtype: {list[DeconflictedRoute]}
        """

        # TODO: for future, should add removal of any time conflicts between the routes (e.g. transition times or temporal keep out zones around the activities in self.routes). Would probably add time keep out zones in for dr_1 in self.routes: loop and check them in for dr_2_indx, dr_2 in enumerate (rr_other.routes): loop. (search the todos in this file!)

        deconflicted_routes = []

        #  if self has routes, then it's possible that these routes share some of the same data from the routes in the other route record ( on the other satellite). this could come up because the routes share some of the same initial windows. we need to determine where the routes split, and use that to determine how much data volume is actually unique ( hasn't already arrived on self), and available to receive
        if len(self._routes) > 0:

            run_optimization =  True

            # figure out window data volume availability per the currently selected routes in self
            avail_dv_by_wind = {}
            for dr_1 in self._routes:
                for wind in dr_1:

                    # if this window is the observation, we want to allow more data volume to be spoken for. This is because it's best to allow more througput than the nominal observation dv to be routed to another other satellites so that if there's contention for a set of windows amongst multiple observations, then there are more routes to provide more choices for pushing data through
                    if type(wind) == ObsWindow:
                        avail_dv_by_wind.setdefault(wind,wind.data_vol * routable_obs_dv_multiplier)
                    else:
                        avail_dv_by_wind.setdefault(wind,wind.data_vol)

                    # subtract off the routes data volume usage from the window's capacity
                    avail_dv_by_wind[wind] -= dr_1.data_vol

            #  record any new windows in the other routes, and for each record which routes use those windows
            other_rts_by_wind = {}
            
            # list of indices that we'll need to drop because of limited data volume availability
            # we mostly just search in this list to see if an element has already been dropped. not making it into a set because I expect the list to be small
            dr_2_indcs_drop = []

            #  now add in the windows found in the other routes, checking as we go if there is enough data volume availability for them to actually be included in the deconflication optimization below
            for dr_2_indx, dr_2 in enumerate (rr_other._routes):
                # todo: to implement transition time constraints, I'd probably add a check right in here that the new xlnk window being appended to the route starts sufficiently far after the current end of the route
                # todo: for overlap, I'd check for any disallowed overlap with any of the winds in any of self._routes

                keep_this_indx = True
                for wind in dr_2:
                    #  if we didn't yet encounter this window in any of the routes in self
                    avail_check = avail_dv_by_wind.setdefault(wind,wind.data_vol)

                    #  this dictionary is used for constraints in the optimization. only want to make a constraint from the window if there's actually sufficient data volume left in the window to select a route through it
                    if avail_check >= min_dv:
                        # record the fact that this window is contained in dr_2
                        other_rts_by_wind.setdefault(wind,[]).append(dr_2_indx)
                    else:
                        #  if there's not enough availability for this window based on the routes in self, add to drop list
                        if not dr_2_indx in dr_2_indcs_drop:
                            dr_2_indcs_drop.append(dr_2_indx)

            #  now that we've sorted out which data routes from rr_other use which Windows, let's make another pass and see if there are data volume conflicts between the routes
            for windex, (wind,dr_2_indcs) in  enumerate (other_rts_by_wind.items()):
                num_rts_possible_wind = floor(avail_dv_by_wind[wind]/min_dv)
                if num_rts_possible_wind < len(dr_2_indcs):
                    # TODO: should some kind of sorting be added here, to give preferential treatment to certain data routes within dr_2_indcs?

                    # mark all data routes pass the maximum number of possible routes through this window for dropping later. 
                    # note the implicit lexicographic preference,  i.e. data routes with earlier indices are dropped last
                    
                    #  update the list of indices based on those that have already been dropped
                    dr_2_indcs_update = [indx for indx in dr_2_indcs if not indx in dr_2_indcs_drop]
                    #  update the drop list if there's anything more that needs to be dropped
                    dr_2_indcs_drop += dr_2_indcs_update[num_rts_possible_wind-1:]


            #  if it turns out no more routes can be selected...
            if len(dr_2_indcs_drop) == len(rr_other._routes):
                run_optimization = False

            

            # construct the data structures for linear program solution
            #  - decision variables are how much data volume to allocate to each route
            #  -- bounds on decision variables, constrained by minimum required route volume
            #  - objective function, cost :  sum of data volumes across all routes ( equally weighted for all routes)
            #  - constraints, A_ub and b_ub:  constrain the sum of the data volumes for all routes going through a given window to be less than or equal to the available data volume for that window

            num_vars = len(rr_other._routes)
            bounds = []
            cost = []
            b_ub = []
            A_ub = []

            #  construct cost and bounds
            # note that costs are equal because we're assuming data volume from any route is the same. If this assumption changes, then this overall algorithm will break (particularly due to rt_possible_cntr below).
            for dr_2_indx in range(num_vars):
                if dr_2_indx in dr_2_indcs_drop:
                    #  if not considering this index, ignore cost and bounds
                    cost.append(0)
                    #  let the data route volume bottom out at zero
                    bounds.append((None,None))
                else:
                    #  negative one because the sense of linprog is minimization
                    cost.append(-1)
                    #  bound on minimum data volume, no upper bound, that's handled by the constraints
                    dr_2_dv = rr_other._routes[dr_2_indx].data_vol
                    bounds.append((min_dv,dr_2_dv))

            #  construct constraints from every window included on selectable routes
            for windex, (wind,dr_2_indcs) in  enumerate (other_rts_by_wind.items()):

                #  we will only keep the window if not all the routes that pass through it have been dropped
                keep_wind = False

                A_ub_row = [0 for dr_2_indx in range(num_vars)]
                #  only factor in those data routes that cross through this window and have not been dropped
                for dr_2_indx in dr_2_indcs:
                    if dr_2_indx in dr_2_indcs_drop: 
                        pass
                    else:
                        A_ub_row[dr_2_indx] = 1
                        keep_wind = True

                if keep_wind:
                    # append the row that specifies which data routes are contributing to the window
                    A_ub.append(A_ub_row)
                    #  constraint imposed by data volume availability from a single window
                    b_ub.append(avail_dv_by_wind[wind])

            #  if all of the windows were dropped because there is no availability, then we can't schedule anything
            if len(b_ub) == 0:
                run_optimization = False


            if run_optimization:
                #  run the solver to figure out the best allocation of data volume to the route possibilities
                res = linprog(cost, A_ub=A_ub, b_ub=b_ub, bounds = bounds,method='simplex')
                
                #  we were not able to find a solution -  either number of iterations ran up, or there is not enough availability to have any de-conflicted routes
                if not res['status'] == 0:
                    # for time being, error out...
                    raise RuntimeWarning('LP not solved successfully')
                    # return []

                #  extract the selected solutions as de-conflicted routes
                dr_2_dv_avails = res['x']
                for dr_2_indx, dr_2 in enumerate (rr_other._routes):
                    dr_dv = dr_2_dv_avails[dr_2_indx]
                    # if dr_2_indx in dr_2_indcs_keep and dr_dv >= min_dv: 
                    if (not dr_2_indx in dr_2_indcs_drop) and dr_dv >= min_dv: 
                        deconflicted_routes.append(DeconflictedRoute(dr=dr_2,available_dv=dr_dv))

            if verbose:
                print('get_dv_deconflicted_routes()')
                if run_optimization:
                    print('  status %d'%(res['status']))
                    print('  num simplex iters %d'%(res['nit']))
                else:
                    print('  did not run optimization')
                print('  num other routes %d'%(len(rr_other._routes)))
                # print('  num routes possible %d'%(len(dr_2_indcs_keep)))
                print('  num routes possible %d'%(num_vars-len(dr_2_indcs_drop)))
                print('  num deconflicted routes %d'%(len(deconflicted_routes)))


        #  if self has no routes ( it hasn't yet received data from anywhere), then every route in  the other route record (on the other satellite) is valid
        else:
            for dr_2 in rr_other._routes:
                available_dv = dr_2.data_vol
                if available_dv >= min_dv: 
                    deconflicted_routes.append(DeconflictedRoute(dr=dr_2,available_dv=available_dv)) 

        return deconflicted_routes


class GPDataRouteSelection():
    """docstring for GP route selection"""

    def __init__(self,gp_params):
        """initializes based on parameters
        
        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """

        scenario_params = gp_params['orbit_prop_params']['scenario_params']
        sat_params = gp_params['orbit_prop_params']['sat_params']
        rs_general_params = gp_params['gp_general_params']['route_selection_general_params']
        rs_v2_params = gp_params['gp_general_params']['route_selection_params_v2']
        as_params = gp_params['gp_general_params']['activity_scheduling_params']
        gp_general_other_params = gp_params['gp_general_params']['other_params']
        link_params = gp_params['orbit_link_params']['general_link_params']
        gp_inst_planning_params = gp_params['gp_instance_params']['planning_params']
        orbit_params = gp_params['orbit_prop_params']['orbit_params']


        self.gp_agent_ID = gp_params['gp_instance_params']['gp_agent_ID']

        # If including cross-links are not when creating routes
        self.include_crosslinks =rs_general_params['include_crosslinks']

        self.num_sats=sat_params['num_sats']
        #  the end of the route selection search window for a given obs will either be the time input from the instance params file, or the end time of the obs plus the filter window length
        self.planning_start_dt  = gp_inst_planning_params['planning_start_dt']
        self.planning_end_obs_dt  = gp_inst_planning_params['planning_end_obs_dt']
        self.planning_end_xlnk_dt  = gp_inst_planning_params['planning_end_xlnk_dt']
        self.planning_end_dlnk_dt  = gp_inst_planning_params['planning_end_dlnk_dt']

        # get the smallest time step used in orbit link. this is the smallest time step we need to worry about in data route selection
        self.act_timestep = min(link_params['xlnk_max_len_s'],link_params['dlnk_max_len_s'])

        #  minimum step1 route data volume -  the minimum data volume a route may have to be considered as a candidate for moving data volume from an observation to a downlink,  when finding routes in step one
        self.min_rs_route_dv =rs_general_params['min_rs_route_dv_Mb']

        # the minimum additional data volume that a data route must be able to add to a data multi route to be considered for adding to that data multi route. Used in accumulate_dr() in routing_objects.py. Not super important, but should be large enough that a DMR doesn't pointlessly add new DRs to it if they don't provide a significant amount of additional DV (considering they add new constraints due to additional windows in the DMR)
        self.min_dmr_candidate_dv = self.min_rs_route_dv / 2  

        # minimum data volume that is considered for routes in the activity scheduling stage
        self.min_obs_dv_dlnk_req =as_params['min_obs_dv_dlnk_req_Mb']

        self.step2_params = rs_v2_params['step2_params']
        
        #  the "effectively zero" number.
        self.dv_epsilon = as_params['dv_epsilon_Mb']

        self.latency_params = gp_params['gp_general_params']['other_params']['latency_calculation']

        self.final_route_records = None

        self.act_timing_helper = ActivityTimingHelper(sat_params['activity_params'],orbit_params['sat_ids_by_orbit_name'],sat_params['sat_id_order'],None) 

        # specifies how much data volume from a given obs is allowed to be selected for routing to other satellites. Want this to be greater than one so that routes account for more than just the exact amount of the obs dv, so that there's more choice in routes to ground 
        self.routable_obs_dv_multiplier = 8

        self.max_num_dlnks_allowed_after_planning_end_xlnk = gp_inst_planning_params['max_num_dlnks_allowed_after_planning_end_xlnk']



    def filter_windows(self,dlnk_winds_flat,xlnk_winds,num_sats,start,dlnk_end_dt,xlnk_end_dt,trim_windows_at_start=False):
        # todo: this should be updated to use gp_general_tools

        dlink_winds_flat_filtered = [[] for sat_indx in  range (num_sats)]
        xlink_winds_flat_filtered = [[[] for xsat_indx in  range ( num_sats)] for sat_indx in  range (num_sats)]
        wind_ids_seen = set()

        for sat_indx in  range (num_sats):
            for xsat_indx in  range ( num_sats):
                for wind in xlnk_winds[sat_indx][xsat_indx]:
                    #  filter out any redundant windows created by a possible earlier copying somewhere
                    if wind.window_ID in wind_ids_seen:
                        continue

                    if wind.duration.total_seconds() < self.act_timing_helper.get_act_min_duration(wind):
                        continue

                    # min_end = min(end_dt_by_sat_indx[sat_indx],end_dt_by_sat_indx[xsat_indx])
                    if  wind.original_start >= start  and  wind.original_end  <= xlnk_end_dt:
                        xlink_winds_flat_filtered[sat_indx][xsat_indx]. append ( wind)
                        wind_ids_seen.add(wind)

            for wind in dlnk_winds_flat[sat_indx]:
                #  filter out any redundant windows created by a possible earlier copying somewhere
                if wind.window_ID in wind_ids_seen:
                    continue

                if wind.duration.total_seconds() < self.act_timing_helper.get_act_min_duration(wind):
                    continue

                if  wind.original_start >= start  and  wind.original_end  <= dlnk_end_dt:
                    dlink_winds_flat_filtered[sat_indx]. append ( wind)
                    wind_ids_seen.add(wind)

                # todo: should clean up this code at some point - it's no longer needed
                # Consider case where the start overlaps with the window, but the center of the window is past the start so we can still get some data volume from the window.  we do this to allow down links that are overlapping with an observation to be considered for that observation -  in practice it turns out to be a large sacrifice to not allow such dumplings to execute ( dictation put dumplings instead of down links, but I'm just gonna leave that there :D WUBBA LUBBA DUB DUB)
                # elif trim_windows_at_start and (wind.original_start < start and wind.center > start and wind.original_end  <= dlnk_end_dt):
                #     wind_copy = deepcopy(wind)
                #     wind_ids_seen.add(wind)
                #     wind_copy.original_wind_ref = wind
                #     wind_copy.modify_time(start,'start')  #  update start and end time. also updates data volume
                #     dlink_winds_flat_filtered[sat_indx]. append ( wind_copy)

                # Consider case where the start overlaps with the window, but the center of the window is past the start so we can still get some data volume from the window.  we do this to allow down links that are overlapping with an observation to be considered for that observation -  in practice it turns out to be a large sacrifice to not allow such dumplings to execute ( dictation put dumplings instead of down links, but I'm just gonna leave that there :D WUBBA LUBBA DUB DUB)
                elif wind.original_start < start and wind.center > start and wind.original_end  <= dlnk_end_dt:
                    wind_copy = deepcopy(wind)
                    wind_ids_seen.add(wind)
                    wind_copy.modify_time(start,'start')  #  update start and end time. also updates data volume
                    dlink_winds_flat_filtered[sat_indx]. append ( wind_copy)

        return dlink_winds_flat_filtered, xlink_winds_flat_filtered

    @staticmethod
    def get_best_rr(rr_dancecards,act,tp_indx,sat_indx,min_route_dv,act_timing_helper):
        """  todo: update this description return the best, valid route record before the start of an activity. The start of an activity ("act") can happen anywhere between two time points. we want to check the second time point to see if there was an activity that ended before act and updated the route record. if not, use the first time point, because any route record on that time point should definitely have ended before act. the use of this function is a means to get around the fact that two cross-links which overlap a single timestep but are temporally consistent ( the second starts after the first ends) are modeled with the first cross-link ending on the time point after this time step, in the second cross-link starting on the time point before this time step"""
        if tp_indx == 0:
            raise RuntimeError('Should not be called with tp_indx = 0')

        rr_candidate = None
        # search backward through the route records for sat_indx to find the one that is the best candidate for act
        while tp_indx >= 0:
            rr_candidate = rr_dancecards[sat_indx][tp_indx]
            tp_indx -= 1

            if rr_candidate is None:
                continue

            # use activity center time here so we know that at least some of the activity is executable after the route record candidate
            # TODO: this right here needs update
            if act.center >= rr_candidate.release_time:
                routable_dv = rr_candidate.get_act_dv_availability(act,sat_indx,act_timing_helper)

                # if we've found a route record that allows us to accept incoming data volume from act, go ahead and return that candidate
                if routable_dv > min_route_dv:
                    return rr_candidate,routable_dv

        return None,0

    def run_step1 ( self,obs_wind,dlnk_winds_flat,xlnk_winds, verbose = False):
        #  note that a potential source of slowness in the code below is the creation of new RouteRecord objects for every sat at every time step 
        # TODO: figure out if there's a more efficient way to do this -  for example, we shouldn't have to preserve these objects when they're more than a certain number of time steps in the past

        # clear state before we run
        self.final_route_records = None

        dr_uid = 0


        # print("ids: dlnk_winds_flat %d xlnk_winds %d"%(id(dlnk_winds_flat),id(xlnk_winds)))

        start_dt = obs_wind.center
        # planning_end_dlnk_dt should be > planning_end_obs,xlnk_dt, to allow the sats to reach into the future for dlnk "backhaul" - high latency, bulk DV delivery
        dlnk_end_dt = self.planning_end_dlnk_dt
        xlnk_end_dt = self.planning_end_xlnk_dt

        # crazy looking line, but it's easy... dictionary of end times by sat_indx - end_dt if not observing sat, else end_obs_sat_dt
        # end_dt_by_sat_indx = {sat_indx: end_dt if sat_indx != obs_wind.sat_indx else end_obs_sat_dt for sat_indx in range (self.num_sats)}
        

        dlnk_winds_flat_filt,xlnk_winds_filt =  self.filter_windows (dlnk_winds_flat,xlnk_winds, self.num_sats, start_dt, dlnk_end_dt, xlnk_end_dt , trim_windows_at_start=True)

        # if obs_wind.window_ID == 43:
        #     debug_tools.debug_breakpt()

        if verbose:
            print ('Running route selection for obs: %s'%(obs_wind))
            print ('from %s to %s'%(start_dt, dlnk_end_dt))


        # construct a set of dance cards for every satellite, 
        # each of which keeps track of all of the activities of satellite 
        # can possibly execute at any given time slice delta T. 
        #  start the dance card from the end of the observation window, because we are only selecting routes after the observation
        act_dancecards = [Dancecard(start_dt,dlnk_end_dt,self.act_timestep,mode='timestep') for sat_indx in range (self.num_sats)]
        
        #  make another set of dance cards that keep track of the route record data structures at a given time point
        #  note that even though these dance cards are in different modes, they share the same time points and time steps (both values and indices)
        rr_dancecards = [Dancecard(start_dt,dlnk_end_dt,self.act_timestep,item_init=None,mode='timepoint') for sat_indx in range (self.num_sats)]

        for sat_indx in range (self.num_sats): 
            act_dancecards[sat_indx].add_winds_to_dancecard(dlnk_winds_flat_filt[sat_indx])

            if self.include_crosslinks:
                #  cross-link Windows matrix is symmetric ( and upper triangular)
                for xsat_indx in  range (sat_indx+1, self.num_sats):
                    act_dancecards[sat_indx].add_winds_to_dancecard(xlnk_winds_filt[sat_indx][xsat_indx])
                    act_dancecards[xsat_indx].add_winds_to_dancecard(xlnk_winds_filt[sat_indx][xsat_indx])
            #  first point for every non-observing satellite has a zero data volume, no path object
            rr_dancecards[sat_indx][0] = RouteRecord(dv=0,release_time = start_dt,routes=[])

        allowed_obs_dv = obs_wind.data_vol * self.routable_obs_dv_multiplier

        #  observing sat starts with an initial data route that includes only the observation window
        #  mark this first data route is having the observation data volume with multiplier included. While there is not actually more observation data volume available to route then is present in the observation window, we can still allow the route to be larger. constraints on actual routed data volume will be enforced at activity scheduling
        first_dr = DataRoute(agent_ID=self.gp_agent_ID,agent_ID_index=dr_uid, route =[obs_wind], window_start_sats={obs_wind: obs_wind.sat_indx},dv=allowed_obs_dv,dv_epsilon =self.dv_epsilon,obs_dv_multiplier=self.routable_obs_dv_multiplier)
        dr_uid += 1
        #  put this initial data route on the first route record for the observing satellite
        rr_dancecards[obs_wind.sat_indx][0] = RouteRecord(dv=allowed_obs_dv,release_time = start_dt,routes=[first_dr])

        #  all dancecards here share the same intial timepoint indices (including for sat_indx, though it's longer)
        time_getter_dc = act_dancecards[obs_wind.sat_indx]

        def time_within(t1,t2,toi):
            return toi >= t1 and toi < t2

        final_route_records = []

        visited_act_set = set()

        #  the seconds value of the very first time point
        tp_dt = time_getter_dc.get_tp_from_tp_indx(0,out_units='datetime')


        num_dlnks_found_after_planning_end_xlnk_by_sat_indx = [0 for sat_indx in range(self.num_sats)]        
        after_planning_end_xlnk = False

        #  the main loop for the dynamic programming algorithm
        # note: have to get the generator again
        for tp_indx in time_getter_dc.get_tp_indcs ():
            if verbose:
                if tp_indx % 10 == 0:
                    print(('tp_indx: %d/%d')%(tp_indx,time_getter_dc.num_timepoints-1))

            #  nothing happens at the first time point index, because that's right at the end of the observation window
            if tp_indx == 0:
                continue


            # get time point values
            tp_last_dt = tp_dt
            tp_dt = time_getter_dc.get_tp_from_tp_indx(tp_indx,out_units='datetime')

            # if after xlnk planning end
            if tp_dt > self.planning_end_xlnk_dt:
                after_planning_end_xlnk = True

            for sat_indx in range (self.num_sats):
                #  if we reach the end of time for this sat index, go to next satellite (i.e. were only looking for down links from observing sat at this point)
                if tp_dt > dlnk_end_dt:
                    continue

                #  get the activities that were active during the time step immediately preceding time point
                acts = act_dancecards[sat_indx].get_objects_at_ts_pre_tp_indx(tp_indx)

                #  get last route record
                rr_last = rr_dancecards[sat_indx][tp_indx-1]

                xlnk_options =[]

                ################
                # go through all the activities active for this time step. Make a list of cross-link windows that are candidates for sending data to sat_indx
                for act in acts:
                    # if we have already considered this activity, keep going ( not currently using for cross-links but leaving here for consistency)
                    # Technically this shouldn't be necessary, because a down/crosslink should only have a start/end time within one time step. but performing this type of check here should be more efficient than checking the start time again
                    if act in visited_act_set: continue

                    if type(act) == XlnkWindow:
                        #  if the cross-link ends in the time step leading up to the current time point, and sat_indx is a receiving satellite in the cross-link
                        if time_within(tp_last_dt,tp_dt,act.original_end) and act.is_rx(sat_indx):
                            xlnk_options.append(act)
                            # visited_act_set.add(act)

                ################
                #  determine for every cross-link option how much data volume the cross-link can deliver to  satellite sat_indx.

                #  each xlnk_candidate is a single xlnk window. For each of these candidates, there's a set of de-conflicted routes that move up to xlnk.data_vol amount of data volume through that xlnk window
                xlnk_candidates = []
                best_dv_seen = 0
                for xlnk in xlnk_options:
                    #  remove redundant calculations: if we've already found a cross-link option that can transport more data volume than  the maximum possible for the current cross-link, skip the current cross-link
                    if xlnk.data_vol <= best_dv_seen:
                        continue



                    #  figure out which route record corresponded to the time right before this cross-link started. we can't assume it's only one time point in the past, because the cross-link could stretch across multiple timesteps
                    tp_indx_pre_xlnk = time_getter_dc.get_tp_indx_pre_t(xlnk.original_start,in_units='datetime')

                    # if xlnk.window_ID == 1500 and obs_wind.window_ID == 39:
                    #     debug_tools.debug_breakpt()

                    # add one to the tp_indx, because we also want to consider any route record with a release time between the timepoints at tp_indx_pre_xlnk and tp_indx_pre_xlnk+1 (because the release could be not-exactly-on-a-tp-indx)
                    rr_last_sat,avail_dv_sat = self.get_best_rr(rr_dancecards,xlnk,tp_indx_pre_xlnk+1,sat_indx,self.min_rs_route_dv,self.act_timing_helper)
                    
                    #  get the route record for the corresponding crosslink partner satellite
                    xsat_indx=xlnk.get_xlnk_partner(sat_indx)
                    rr_last_xsat,avail_dv_xsat = self.get_best_rr(rr_dancecards,xlnk,tp_indx_pre_xlnk+1,xsat_indx,self.min_rs_route_dv,self.act_timing_helper)


                    avail_dv_for_xlnk = min(avail_dv_sat,avail_dv_xsat)


                    if rr_last_sat is None or rr_last_xsat is None or avail_dv_for_xlnk < self.min_rs_route_dv:
                        # this shouldn't occur very often at all - the only case I can think of is if rr_last_xsat is right after the observation on the obs sat, and xlnk overlaps with that obs
                        continue

                    #  remove redundant calculations: if the other satellite has less data volume at this time point than we do, then we can't get any more data volume from them, and it's pointless to consider it
                    if rr_last_sat.dv > rr_last_xsat.dv:
                        continue

                    #  need to figure out what data on the other satellite is data that we have not yet received on sat_indx. this returns a set of de-conflicted routes that all send valid, non-duplicated data to sat_indx
                    # note that this deconfliction does not currently look at transition times across routes, and so has some slop...
                    deconf_rts = rr_last_sat.get_dv_deconflicted_routes(rr_last_xsat,min_dv=self.min_rs_route_dv,routable_obs_dv_multiplier=self.routable_obs_dv_multiplier,verbose= False)

                    #todo: hmm, do we also want to consider other direction, where rr_sat is incumbent and rr_last_sat are the routes to deconflict?

                    #  now we need to figure out how many of these routes we can use, based upon the available crosslink bandwidth
                    xlnk_candidate_rts = []
                    # Cumulative data volume used on the crosslink by the route candidates. must be less than the available crosslink data volume
                    x_cum_dv = 0
                    # deconf_rt is a DeconflictedRoute namedtuple, from above
                    for deconf_rt in deconf_rts:
                        # if we're less than the total crosslink data volume
                        if x_cum_dv + deconf_rt.available_dv <= avail_dv_for_xlnk:
                            xlnk_candidate_rts.append(deconf_rt)
                            x_cum_dv += deconf_rt.available_dv
                        # also handle the case where we don't have enough available data volume to fulfill all of the potential for route, but we can still grab enough data volume for that route ( greater than specified minimum) to use it
                        else:
                            available_dv = avail_dv_for_xlnk - x_cum_dv 
                            #  check if greater than specified minimum
                            if available_dv >= self.min_rs_route_dv:
                                deconf_rt.available_dv = available_dv
                                xlnk_candidate_rts.append(deconf_rt)
                                x_cum_dv += available_dv


                    #  make a new candidate entry with a record of the new data volume that we'll have if we choose that candidate.  also bring along some other relevant objects for bookkeeping
                    new_dv = rr_last_sat.dv + x_cum_dv
                    best_dv_seen = max(best_dv_seen,x_cum_dv)
                    # XLNK_CANDIDATES:       0         1            2                3
                    xlnk_candidates.append((new_dv,rr_last_sat,xlnk_candidate_rts,xlnk))


                found_candidates = len(xlnk_candidates) > 0

                ################
                #  here is where we implement the Bellman equation.
                #  find the cross-link that delivers the most data volume
                if found_candidates:
                    # 0 indx is "new_dv"
                    best_xlnk_cand = max(xlnk_candidates,key= lambda cand: cand[0])


                ################
                #  now, for the best cross-link candidate (if there is one),  make a new route record
                if found_candidates and best_xlnk_cand[0] > rr_last.dv:

                                    
                    # Note: the indices called out below correspond to the fields labeled in XLNK_CANDIDATES line above
                    #  need to make a copy here because the last route record may be visited again, and the object needs to be left as is
                    rr_new = copy(best_xlnk_cand[1])
                    rr_new.dv = best_xlnk_cand[0]
                    xlnk_candidate_rts = best_xlnk_cand[2]
                    xlnk_wind = best_xlnk_cand[3]
                    #  the release time of this new route record is the end of the cross-link - i.e. this route record is valid for all t past the end of the cross-link
                    #  note that it may be possible that the "best" candidate cross-link (most dv moved) has a later release time than another candidate, and this causes some later cross-links to be ruled out because they start before the best (and end after the other candidate). We'll assume this is an acceptable error though. the time step for the dance cards should be made small enough that this is not a big deal.
                    rr_new.release_time = xlnk_wind.center

                    # deconf_rt is a DeconflictedRoute namedtuple, from above
                    for deconf_rt in xlnk_candidate_rts:
                        #  again, make a copy because the existing data routes need to be left as is for future consideration in the algorithm
                        new_dr = copy(deconf_rt.dr)
                        #  the available data volume is simply that of the de-conflicted route, because we already made sure above that the sum of the available data volume from all de-conflicted routes is less than or equal to the throughput of the cross-link, and none of that data volume conflicts  with data volume already present in rr_new
                        new_dr.data_vol = deconf_rt.available_dv
                        new_dr.set_id(self.gp_agent_ID,dr_uid)
                        #  we did not append the cross-link window to the route before ( for sake of efficiency), so now we need to do it. The second argument below records the fact that the cross-link started on the cross-link partner,  not on satellite sat_indx
                        new_dr.append_wind_to_route(xlnk_wind,window_start_sat_indx=xlnk_wind.get_xlnk_partner(sat_indx))
                        
                        rr_new.add_to_routes(new_dr)
                        dr_uid +=1
                else:
                    #  need to make a copy here because the last route record may be visited again, and needs to be left as is
                    rr_new = copy(rr_last)
                
                #  update the dance card at current time point.
                # this makes a new route record at the time point after the time step in which a cross-link ended. It is not necessarily the case that a crosslink ended exactly on a timepoint
                rr_dancecards[sat_indx][tp_indx] = rr_new


                ################
                # go through all the activities active for this time step. pinch off any routes that can end at downlinks.
                #  wanted to wait until after performing all cross-links to check this, because cross-links may have delivered some data volume before the downlink starts in the midst of this timestep
                for act in acts:

                    # if we have already considered this activity, keep going
                    # Technically this shouldn't be necessary, because a down/crosslink should only have a start/end time within one time step. but performing this type of check here should be more efficient than checking the start time again
                    if act in visited_act_set: continue

                    if type(act) == DlnkWindow:


                        # don't consider a dlnk if its center time doesn't fall within or after this timestep (assumption: center time of the activity must be within the final time window of execution)
                        if not act.center > tp_last_dt:
                            continue

                        # check if we've already reached our quota for num allowed dlnks after the xlnks (if enforced)
                        if self.max_num_dlnks_allowed_after_planning_end_xlnk is not None:
                            if (after_planning_end_xlnk and 
                                num_dlnks_found_after_planning_end_xlnk_by_sat_indx[sat_indx] > self.max_num_dlnks_allowed_after_planning_end_xlnk):
                                continue

                        # we have found a dlnk in this timestep, but in order to use it, we need it to START in this timestep. So if it doesn't already, update its start time so it will. We need to make a deepcopy of the window so this change doesn't step on the toes of other routes
                        dlnk = act # save before copying
                        # deal with case where downlink actually started before this timestep
                        if act.start < tp_last_dt:
                            dlnk = deepcopy(act)
                            new_start = tp_last_dt 
                            dlnk.modify_time(new_start,'start')  #  update start and end time. also updates data volume

                        #  figure out if we want to use the route record from the last timepoint, or if a cross-link has delivered more data volume. grab the latter if valid
                        rr_pre_dlnk,avail_dv_for_dlnk = self.get_best_rr(rr_dancecards,dlnk,tp_indx,sat_indx,self.min_rs_route_dv,self.act_timing_helper)

                        # if obs_wind.window_ID == 43:
                        #     if dlnk.window_ID == 81:
                        #         debug_tools.debug_breakpt()

                        # couldn't find any routes valid for this dlnk to send down
                        if rr_pre_dlnk is None or avail_dv_for_dlnk < self.min_rs_route_dv:
                            continue

                        # todo: should check transition time between last xlnk in each of the routes and the start of the downlink - I'd add that check here, or maybe in the loop below

                        #  if the route record actually shows some data volume arriving at this satellite, AND we have a downlink, then we have found an optimal route to the downlink. Create a final route record and save for later.
                        if rr_pre_dlnk.dv > 0:
                            # rr_dlnk = copy(rr_pre_dlnk)
                            #  note: release time no longer matters because were not putting this route record back into the dance card.
                            rr_dlnk = RouteRecord(dv=rr_pre_dlnk.dv,release_time=dlnk.center,routes=[])
                            
                            # available data volume is limited by how much we had at the last route record ( which could be as much as the observation data volume multiplied by the routable_obs_dv_multiplier factor)
                            available_dv = rr_pre_dlnk.dv

                            #  now we need to update the routes within the route record to include the downlink window
                            for dr in rr_pre_dlnk:
                                #  for each route, we will only allocate as much data volume as is deliverable -  the minimum of:
                                # - the data volume that we found for the route ( which is constrained by all of the cross-link capacities in the route)
                                # - the available data volume from this route record
                                # - the bservation data volume
                                # - the downlink data volume
                                # note that it is not necessary in activity scheduling for a data route to conform to all of these limits, because capacity constraints are checked for all of the windows within all the routes. nonetheless, it's good practice to only allocate as much data volume as is really present to every route
                                dv_slice = min(dr.data_vol,available_dv,dr.get_obs().data_vol,avail_dv_for_dlnk)

                                new_dr = copy(dr)
                                new_dr.data_vol = dv_slice
                                new_dr.set_id(self.gp_agent_ID,dr_uid)
                                new_dr.append_wind_to_route(dlnk,window_start_sat_indx=sat_indx)
                                new_dr.fix_window_copies()
                                rr_dlnk.add_to_routes(new_dr)

                                dr_uid +=1
                                available_dv -= dv_slice

                                #  we run out of available data volume so no more routes can be grabbed
                                if available_dv < self.dv_epsilon:
                                    break

                            #  sanity check:  make sure that data volume for the route record greater than or equal to the sum of the data volumes for all of the routes
                            assert(rr_dlnk.dv + self.dv_epsilon >= sum(dr.data_vol for dr in rr_dlnk))

                            rr_dlnk.dv = sum(dr.data_vol for dr in rr_dlnk)

                            final_route_records.append(rr_dlnk)
                            visited_act_set.add(dlnk)

                            # mark as being after end of xlnk planning time, if applicable
                            if after_planning_end_xlnk:
                                num_dlnks_found_after_planning_end_xlnk_by_sat_indx[sat_indx] += 1


        self.final_route_records = final_route_records

        all_routes = [dr for rr in final_route_records for dr in rr]
        for dr in all_routes:
            dr.validate(self.act_timing_helper,time_option='center')

        return all_routes

    def best_route_dv_availability(self,rts,dv_avail_by_wind,check_availability):
        # if we're not paying attention to which windows have already had dv "spoken for", then just return the first index
        if not check_availability:
            if len(rts) > 0:
                best_rt = rts[0]
                best_rt_indx = 0
                return best_rt,best_rt_indx

        best_rt_dv = None
        best_rt = None
        best_rt_indx = None
        for rt_indx,rt in enumerate(rts):
            rt_dv = rt.data_vol
            for act in rt.get_winds():
                # ignore the data availability of the obs windows
                if type(act) == ObsWindow:
                    continue
                dv_avail_by_wind.setdefault(act,act.data_vol)
                rt_dv = min(rt_dv,dv_avail_by_wind[act])

            if (best_rt_dv and rt_dv > best_rt_dv) or best_rt_dv is None:
                best_rt_dv = rt_dv
                best_rt = rt
                best_rt_indx = rt_indx

        return best_rt,best_rt_indx

    def get_from_sorted(self,sorted_rts,num_rts,min_dmr_candidate_dv,drs_taken,dv_avail_by_wind,existing_routes_set=set(),check_availability=False):
        rt_index = 0
        sel_rts = []
        rts_remaining = sorted_rts
        while len(rts_remaining) > 0:            
            curr_dr,_ = self.best_route_dv_availability(rts_remaining,dv_avail_by_wind,check_availability)

            # don't consider the route if it's already spoken for within another dmr.
            if curr_dr in drs_taken:
                rts_remaining.remove(curr_dr)
                continue

            # If the current dr out we're looking at is an existing route, then it should be already a DataMultiRoute
            if curr_dr in existing_routes_set:
                # ensure type correctness, and that it has sufficient data volume
                assert(type(curr_dr)==DataMultiRoute)
                # WG: only do the assert below if it is a GP made DMR, ignore if made by LP
                if not curr_dr.ID.creator_agent_ID[0] == 'S':
                    assert(curr_dr.data_vol >= self.min_obs_dv_dlnk_req - self.dv_epsilon)
                dmr = curr_dr
            #  create a new data multi-route encapsulating the data route. use the existing route's ID -  it won't need it anymore.
            else:    
                dmr = DataMultiRoute(curr_dr.ID,data_routes=[curr_dr],dv_epsilon=self.dv_epsilon)

            rts_remaining.remove(curr_dr) 
            for act in curr_dr.get_winds():
                dv_avail_by_wind.setdefault(act,act.data_vol)
                dv_avail_by_wind[act] -= curr_dr.data_vol

            #  if the data route already has enough data volume to meet the minimum requirement for the activity scheduling stage, add it as a selected route
            if dmr.data_vol >= self.min_obs_dv_dlnk_req:
                sel_rts.append(dmr)
                for dr in dmr.data_routes: drs_taken.add(dr)

            #  otherwise, we need to smoosh routes together in order to meet that requirement
            else:
                #  consider additional data routes in the list. add them to the data multi-route to see if we can make minimum data volume requirement
                
                while len(rts_remaining) > 0:
                    next_dr,_ = self.best_route_dv_availability(rts_remaining,dv_avail_by_wind,check_availability)

                    if next_dr in drs_taken or next_dr in existing_routes_set:
                        rts_remaining.remove(next_dr)
                        continue

                    #  one more check to make sure that we're not operating on an existing route. it would be pretty bad to modify an existing route, below
                    assert(dmr not in existing_routes_set)

                    #  add the data route to the data multi-route if possible
                    dmr.accumulate_dr( next_dr,min_dmr_candidate_dv)
                    rts_remaining.remove(next_dr)
                    for act in next_dr.get_winds():
                        dv_avail_by_wind.setdefault(act,act.data_vol)
                        dv_avail_by_wind[act] -= next_dr.data_vol

                    #  check if we are meeting the minimum data volume requirement after the accumulation. if yes include the selected route and  and break out of the loop
                    if dmr.data_vol >= self.min_obs_dv_dlnk_req:
                        sel_rts.append(dmr)
                        for dr in dmr.data_routes: drs_taken.add(dr)
                        break


            if len(sel_rts) >= num_rts:
                break 


        return sel_rts,drs_taken,rts_remaining

    def run_step2(self,routes_by_obs,overlap_cnt_by_route,existing_routes_set=set()):

        rts_by_obs_sorted_overlap = {}
        rts_by_obs_sorted_dv = {}
        rts_by_obs_sorted_lat = {}

        def  latency_getter(dr):
            return dr.get_latency(
                    'minutes',
                    obs_option = self.latency_params['obs'], 
                    dlnk_option = self.latency_params['dlnk']
                )

        for obs,rts in routes_by_obs.items():
            #  sort the routes for each observation by data volume. Prefer routes with more data volume
            rts_by_obs_sorted_dv[obs] = sorted(rts,key=lambda dr: dr.data_vol,reverse=True)
            #  sort the routes for each observation by latency. Prefer routes with lower latency
            rts_by_obs_sorted_lat[obs] = sorted(rts,key=lambda dr: latency_getter(dr))
            #  sort the routes for each observation by overlap count. Prefer routes with less overlap
            rts_by_obs_sorted_overlap[obs] = sorted(rts,key=lambda dr: overlap_cnt_by_route[dr])

            
        # dmr_uid = 0

        # selected data multi routes
        selected_dmrs_by_obs = {}
        drs_taken = set()
        dv_avail_by_wind = {}
        obs_last_to_first = sorted(routes_by_obs.keys(),key=lambda act:act.center,reverse=True)
        for obs in obs_last_to_first:
            selected_dmrs_by_obs[obs] = []

            sel_rts,drs_taken,_ = self.get_from_sorted(rts_by_obs_sorted_dv[obs],self.step2_params['num_rts_sel_per_obs_dv'],self.min_dmr_candidate_dv,drs_taken,dv_avail_by_wind,existing_routes_set,check_availability=False)
            selected_dmrs_by_obs[obs] += sel_rts

            sel_rts,drs_taken,_ = self.get_from_sorted(rts_by_obs_sorted_lat[obs],self.step2_params['num_rts_sel_per_obs_lat'],self.min_dmr_candidate_dv,drs_taken,dv_avail_by_wind,existing_routes_set,check_availability=False)
            selected_dmrs_by_obs[obs] += sel_rts

        curr_num_rts = 0
        num_overlap_sel_rts = self.step2_params['num_rts_sel_per_obs_overlap']
        while curr_num_rts < num_overlap_sel_rts:
            for obs in obs_last_to_first:
                if len(rts_by_obs_sorted_overlap[obs]) == 0:
                    continue

                # note: Don't need to worry about a data route appearing multiple times across the selected routes returned from get_from_sorted, because each data route may be only used once (enforced with drs_taken)
                sel_rts,drs_taken,rts_by_obs_sorted_overlap[obs] = self.get_from_sorted(rts_by_obs_sorted_overlap[obs],1,self.min_dmr_candidate_dv,drs_taken,dv_avail_by_wind,existing_routes_set,check_availability=True)
                selected_dmrs_by_obs[obs] += sel_rts

            curr_num_rts += 1

        return selected_dmrs_by_obs

    def get_stats(self,final_route_records,verbose=False):
        stats = {}
        # stats['num_dlnk_windows'] = sum([len (self.dlnk_winds_flat[sat_indx]) for sat_indx in range (self.num_sats)])
        # stats['num_xlnk_windows'] = sum([len (self.xlnk_winds[sat_indx][xsat_indx]) for xsat_indx in range (self.num_sats) for sat_indx in range (self.num_sats)])

        if verbose:
            for rr_indx, rr in  enumerate (self.final_route_records):
                print ('route record %d, %d routes, sum rts/ dlnk dv: %.0f/%.0f Mb (%.2f%%)'%(rr_indx,len(rr),sum(dr.data_vol for dr in rr),rr._routes[0].get_dlnk().data_vol,100*sum(dr.data_vol for dr in rr)/rr._routes[0].get_dlnk().data_vol))
                for dr_indx, dr in  enumerate (rr):
                    pass
                    # print('  dr %d %s'%(dr_indx,dr))
                    print('  dr %d %s'%(dr_indx,dr.get_dlnk ()))

            # print ( "Obs dv: %f" % ( self.obs_wind.data_vol))
            # print ( "Considering %d downlink windows" % (stats['num_dlnk_windows']))
            # print ( "Considering %d crosslink windows" % (stats['num_xlnk_windows']))

        return stats




