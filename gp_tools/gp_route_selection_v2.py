# Algorithm for creating a set of data routes for a given observation, later fed into the activity scheduling stage
# 
# This version uses a dynamic programming algorithm to optimize selected paths
# 
# @author Kit Kennedy
#
#  note that a path is the same as a route. 

from  datetime import timedelta
from copy import copy
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
from .routing_objects import DataRoute
from .schedule_objects import Dancecard
from .custom_activity_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow

DATE_STRING_FORMAT = 'short'
# DATE_STRING_FORMAT = 'iso'

def short_date_string(dt):
    return dt.strftime("%H:%M:%S")

def date_string(dt):
    if DATE_STRING_FORMAT == 'iso':
        return dt.isoformat()
    if DATE_STRING_FORMAT == 'short':
        return  short_date_string(dt)

class DeconflictedRoute():
    """docstring for RouteRecord"""
    def __init__(self,dr,available_dv):
        self.dr = dr
        self.available_dv = available_dv

    def __repr__( self):
        return 'DeconflictedRoute(%s,%s)'%( self.dr, self.available_dv)

def no_route_conflict(dr1,dr2):
    """ check for route conflicts. implementation specific check to see if a cross-link window in one route both overlaps in time and includes the same satellites as a cross-link window in another route"""

     # go through all of the windows in route 1 and for each one see if it conflicts with a window in route 2
    for wind1 in dr1:
        if not type(wind1) == XlnkWindow:
            continue

        for wind2 in dr2:
            if not type(wind2) == XlnkWindow:
                continue

            #  check for overlap
            #  if before the start, continue to next window in route 2
            if wind1.end <= wind2.start:
                continue
            #  if after the end, we've passed a point where wind1 could overlap with a window in route 2.
            elif wind1.start >= wind2.end:
                break
            # if there's any overlap at all...
            else:
                #  I'm sure there's a better way to express this check
                wind2_indcs = [wind2.sat_indx,wind2.xsat_indx]
                if wind1.sat_indx in wind2_indcs or wind1.xsat_indx in wind2_indcs:
                    return False

    return True

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
        self.routes = routes
        self.release_time = release_time

    def __repr__( self):
        return 'RouteRecord(%s,%s,%s)'%( self.dv, date_string(self.release_time),self.routes)

    def get_deconflicted_routes(self,rr_other,min_dv=0,verbose = False):
        """ determines which routes from other route record can be extended to self route record
        
        determines from the set of routes contained in self which routes in rr_other are able to be extended to the satellite and time point for the self route record. determines this by figuring out what data volume is still available in windows to run routes through, and optimally selecting which routes from rr_other are allocated that data volume. optimal here means the allocation that maximizes the data volume available in the de-conflicted routes returned
        :param rr_other:  the other route record, containing a set of data routes possible to extend to self
        :type rr_other: RouteRecord
        :param min_dv:  minimum data volume allowable for a route to be selectable, defaults to 0
        :type min_dv: number, optional
        :param verbose: output verbosity, defaults to False
        :type verbose: bool, optional
        :returns:  a set of de-conflicted routes with data volume allocations that, if all selected together, would violate no constraints on window data volume capacities across all routes in self
        :rtype: {list[DeconflictedRoute]}
        """

        deconflicted_routes = []

        #  if self has routes, then it's possible that these routes share some of the same data from the routes in the other route record ( on the other satellite). this could come up because the routes share some of the same initial windows. we need to determine where the routes split, and use that to determine how much data volume is actually unique ( hasn't already arrived on self), and available to receive
        if len(self.routes) > 0:

            # figure out window data volume availability per the currently selected routes in self
            avail_dv_by_wind = {}
            for dr_1 in self.routes:
                for wind in dr_1:
                    # subtract off the routes data volume usage from the window's capacity
                    avail_dv_by_wind.setdefault(wind,wind.data_vol)
                    avail_dv_by_wind[wind] -= dr_1.data_vol

            #  record any new windows in the other routes, and record which routes use which windows
            other_rts_by_wind = {}
            dr_2_indcs_keep = []
            for dr_2_indx, dr_2 in enumerate (rr_other.routes):
                keep_this_indx = True
                for wind in dr_2:
                    #  if we didn't yet encounter this window in any of the routes in self
                    avail_check = avail_dv_by_wind.setdefault(wind,wind.data_vol)

                    #  this dictionary is used for constraints in the optimization. only want to make a constraint from the window if there's actually sufficient data volume left in the window to select a route through it
                    if avail_check >= min_dv:
                        # record the fact that this window is contained in dr_2
                        other_rts_by_wind.setdefault(wind,[]).append(dr_2_indx)
                    else:
                        keep_this_indx = False

                #  if every window in the row has sufficient data volume ( to meet minimum requirement) this route is possible to select
                if keep_this_indx:
                    dr_2_indcs_keep.append(dr_2_indx)

            #  if it turns out no more routes can be selected...
            if len(dr_2_indcs_keep) == 0:
                return []

            # construct the data structures for linear program solution
            #  - decision variables are how much data volume to allocate to each route
            #  -- bounds on decision variables, constrained by minimum required route volume
            #  - objective function, cost :  sum of data volumes across all routes ( equally weighted for all routes)
            #  - constraints, A_ub and b_ub:  constrain the sum of the data volumes for all routes going through a given window to be less than or equal to the available data volume for that window
            
            b_ub = [0 for i in range(len(other_rts_by_wind.keys()))]
            A_ub = [[0 for j in range(len (rr_other.routes))] for i in range(len (other_rts_by_wind.keys()))]

            #  construct cost
            # note that costs are equal because we're assuming data volume from any route is the same. If this assumption changes, then this overall algorithm will break (particularly due to rt_possible_cntr below).
            cost = []
            for dr_2_indx in range(len(rr_other.routes)):
                if dr_2_indx in dr_2_indcs_keep:
                    #  negative one because the sense of linprog is minimization
                    cost.append(-1)
                else:
                    #  if not considering this index, ignore cost and bounds
                    cost.append(0)

            #  construct constraints from every window included on selectable routes
            for windex, (wind,dr_2_indcs) in  enumerate (other_rts_by_wind.items()):

                #  constraint imposed by data volume availability from a single window
                b_ub[windex] = avail_dv_by_wind[wind]

                #  only factor in those data routes that cross through this window
                #  Note:  a row can be all zeros if a dr_indx was  disqualified from dr_2_indcs_keep  after one of its windows was already added to other_rts_by_wind.  not great, but has no effect on the LP simplex algorithm
                for dr_2_indx in dr_2_indcs:
                    if dr_2_indx in dr_2_indcs_keep: 
                        A_ub[windex][dr_2_indx] = 1

            #  there is going to be a most constrained window that can only allow a certain number of routes through. this constrains the number of routes with min_dv that we can get from our linear program solution. figure out what this number is
            #  note that we might've been able to figure this out in the code above -  adding it here is mildly hacky.
            # consider b_ub = [ 339.07, 201.58, 111.23], A_ub = [[1,1,1,1],[0,1,0,1],[1,0,0,0]], and min_dv = 100. In this case, num_routes_possible is 3, from the 339.07 index, because that most constrains the DESIRED # of routes through a window
            # zero_protect where a row of the A matrix is all zeros.  see note above about this. 
            zero_protect = 0.001
            if len (rr_other.routes) > 1/zero_protect:
                raise RuntimeWarning('zero_protect is not small enough')
            num_routes_possible_indx = [floor(b_ub[indx]/min_dv) for indx in range(len(b_ub))]
            num_routes_desired_indx = [sum(A_ub[indx])+zero_protect for indx in range(len(b_ub))]
            mindex = argmin([num_routes_possible_indx[indx]/num_routes_desired_indx[indx] for indx in range(len(b_ub))])
            num_routes_possible = num_routes_possible_indx[mindex]

            #  construct bounds
            bounds = []
            rt_possible_cntr = 0
            for dr_2_indx in range(len(rr_other.routes)):
                # if we're keeping the data route as a possibility because both:
                # - the dr has no window along its route that has less available data volume than minimum dv
                # - the counter for the number of possible data routes has not reached its limit
                # note that we're implicitly deciding to choose the first routes we see as those that will get selected. This is only okay because all the costs are equal 
                if dr_2_indx in dr_2_indcs_keep and rt_possible_cntr < num_routes_possible:
                    #  bound on minimum data volume, no upper bound, that's handled by the constraints
                    bounds.append((min_dv,None))
                    rt_possible_cntr += 1
                else:
                    #  otherwise, let the data route volume bottom out at zero if need be
                    bounds.append((None,None))

            #  run the solver to figure out the best allocation of data volume to the route possibilities
            res = linprog(cost, A_ub=A_ub, b_ub=b_ub, bounds = bounds,method='simplex')

            if verbose:
                print('get_deconflicted_routes()')
                print('  status %d'%(res['status']))
                print('  num simplex iters %d'%(res['nit']))
                print('  num other routes %d'%(len(rr_other.routes)))
                print('  num routes possible %d'%(len(dr_2_indcs_keep)))

                # if res['status'] == 2:
                #     import ipdb
                #     ipdb.set_trace()

            #  we were not able to find a solution -  either number of iterations ran up, or there is not enough availability to have any de-conflicted routes
            if not res['status'] == 0:
                # for time being, error out...
                raise RuntimeWarning('LP not solved successfully')
                # return []

            #  extract the selected solutions as de-conflicted routes
            dr_2_dv_avails = res['x']
            for dr_2_indx, dr_2 in enumerate (rr_other.routes):
                dr_dv = dr_2_dv_avails[dr_2_indx]
                if dr_2_indx in dr_2_indcs_keep and dr_dv >= min_dv: 
                    deconflicted_routes.append(DeconflictedRoute(dr=dr_2,available_dv=dr_dv))

        #  if self has no routes ( it hasn't yet received data from anywhere), then every route in  the other route record (on the other satellite) is valid
        else:
            for dr_2 in rr_other.routes:
                available_dv = dr_2.data_vol
                if available_dv >= min_dv: 
                    deconflicted_routes.append(DeconflictedRoute(dr=dr_2,available_dv=available_dv)) 

        return deconflicted_routes

    def __copy__(self):
        newone = type(self)( self.dv,copy(self.release_time))
        #  make a shallow copy of every data route object - this avoids copying the underlying windows image data route
        newone.routes = [copy(rt) for rt in self.routes]
        return newone

class GPDataRouteSelection():
    """docstring for GP route selection"""

    def __init__(self,gp_params):
        """initializes based on parameters
        
        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """

        scenario_params = gp_params['gp_orbit_prop_params']['scenario_params']
        sat_params = gp_params['gp_orbit_prop_params']['sat_params']
        rs_general_params = gp_params['gp_general_params']['route_selection_general_params']
        rs_params = gp_params['gp_general_params']['route_selection_params_v2']
        gp_general_other_params = gp_params['gp_general_params']['other_params']
        link_params = gp_params['gp_orbit_link_params']['link_params']
        gp_as_inst_params = gp_params['gp_instance_params']['activity_scheduling_params']

        self.num_sats=sat_params['num_sats']
        #  these times indicate the (largest) window over which we are considering routes
        self.sel_latest_end_dt  = tt.iso_string_to_dt (gp_as_inst_params['end_utc'])

        # get the smallest time step used in orbit link. this is the smallest time step we need to worry about in data route selection
        self.act_timestep = min(link_params['xlnk_max_len_s'],link_params['dlnk_max_len_s'])

        #  minimum simple route data volume -  the minimum data volume a route may have to be considered as a candidate for moving data volume from an observation to a downlink
        self.min_s_route_dv =rs_general_params['min_simple_route_dv_Mb']
        
        #  the "effectively zero" number. defining in terms of the minimum path data volume for future proofing
        self.dv_epsilon = self.min_s_route_dv / 1000
        self.wind_filter_duration =  timedelta (seconds =rs_general_params['wind_filter_duration_s'])
        self.wind_filter_duration_obs_sat =  timedelta (seconds =rs_general_params['wind_filter_duration_obs_sat_s'])
        self.latency_params =  gp_general_other_params['latency_calculation']

        if self.latency_params['obs'] not in ['start','end']:
            raise NotImplementedError
        if self.latency_params['dlnk'] not in ['start','end','center']:
            raise NotImplementedError


    @staticmethod
    def  filter_windows(dlnk_winds_flat,xlnk_winds,num_sats,start,end_dt_by_sat_indx):

        dlink_winds_flat_filtered = [[] for sat_indx in  range (num_sats)]
        xlink_winds_flat_filtered = [[[] for xsat_indx in  range ( num_sats)] for sat_indx in  range (num_sats)]

        for sat_indx in  range (num_sats):
            for xsat_indx in  range ( num_sats):
                for wind in xlnk_winds[sat_indx][xsat_indx]:
                    min_end = min(end_dt_by_sat_indx[sat_indx],end_dt_by_sat_indx[xsat_indx])
                    if  wind.start > start  and  wind.end  < min_end:
                        xlink_winds_flat_filtered[sat_indx][xsat_indx]. append ( wind)

            for wind in dlnk_winds_flat[sat_indx]:
                if  wind.start > start  and  wind.end  < end_dt_by_sat_indx[sat_indx]:
                    dlink_winds_flat_filtered[sat_indx]. append ( wind)

        return dlink_winds_flat_filtered, xlink_winds_flat_filtered

    def run_stage1 ( self,obs_wind,dlnk_winds_flat,xlnk_winds, dr_uid, verbose = False):
        #  note that a potential source of slowness in the code below is the creation of new RouteRecord objects for every sat at every time step 
        # TODO: figure out if there's a more efficient way to do this -  for example, we shouldn't have to preserve these objects when they're more than a certain number of time steps in the past


        start_dt = obs_wind.end
        end_dt = min(self.sel_latest_end_dt,start_dt + self.wind_filter_duration)
        end_obs_sat_dt = min(self.sel_latest_end_dt,start_dt + self.wind_filter_duration_obs_sat)

        # crazy looking line, but it's easy... dictionary of end times by sat_indx - end_dt if not observing sat, else end_obs_sat_dt
        end_dt_by_sat_indx = {sat_indx: end_dt if sat_indx != obs_wind.sat_indx else end_obs_sat_dt for sat_indx in range (self.num_sats)}
        
        dlnk_winds_flat_filt,xlnk_winds_filt =  self.filter_windows (dlnk_winds_flat,xlnk_winds, self.num_sats, obs_wind.end, end_dt_by_sat_indx )

        print ('Running route selection for obs: %s'%(obs_wind))
        print ('from %s to %s'%(start_dt, end_dt))


        # construct a set of dance cards for every satellite, 
        # each of which keeps track of all of the activities of satellite 
        # can possibly execute at any given time slice delta T. 
        #  start the dance card from the end of the observation window, because we are only selecting routes after the observation
        act_dancecards = [Dancecard(start_dt,end_dt_by_sat_indx[sat_indx],self.act_timestep,mode='timestep') for sat_indx in range (self.num_sats)]
        
        #  make another set of dance cards that keep track of the route record data structures at a given time point
        #  note that even though these dance cards are in different modes, they share the same time points and time steps (both values and indices)
        rr_dancecards = [Dancecard(start_dt,end_dt_by_sat_indx[sat_indx],self.act_timestep,item_init=None,mode='timepoint') for sat_indx in range (self.num_sats)]

        obs_dv = obs_wind.data_vol

        for sat_indx in range (self.num_sats): 
            act_dancecards[sat_indx].add_winds_to_dancecard(dlnk_winds_flat_filt[sat_indx])
            #  cross-link Windows matrix is symmetric ( and upper triangular)
            for xsat_indx in  range (sat_indx+1, self.num_sats):
                act_dancecards[sat_indx].add_winds_to_dancecard(xlnk_winds_filt[sat_indx][xsat_indx])
                act_dancecards[xsat_indx].add_winds_to_dancecard(xlnk_winds_filt[sat_indx][xsat_indx])
            #  first point for every non-observing satellite has a zero data volume, no path object
            rr_dancecards[sat_indx][0] = RouteRecord(dv=0,release_time = start_dt,routes=[])

        #  observing sat starts with an initial data route that includes only the observation window
        first_dr = DataRoute(dr_uid, route =[obs_wind], window_start_sats={obs_wind: obs_wind.sat_indx},dv=obs_dv)
        dr_uid += 1
        #  put this initial data route on the first route record for the observing satellite
        rr_dancecards[obs_wind.sat_indx][0] = RouteRecord(dv=obs_dv,release_time = start_dt,routes=[first_dr])

        #  all dancecards here share the same intial timepoint indices (including for sat_indx, though it's longer)
        time_getter_dc = act_dancecards[obs_wind.sat_indx]

        def time_within(t1,t2,toi):
            return toi >= t1 and toi < t2

        def get_best_rr(t,tp_indx_pre_t,sat_indx):
            """  return the best, valid route record before the start of an activity. The start of an activity ("act") can happen anywhere between two time points. we want to check the second time point to see if there was an activity that ended before act and updated the route record. if not, use the first time point, because any route record on that time point should definitely have ended before act. the use of this function is a means to get around the fact that two cross-links which overlap a single timestep but are temporarily consistent ( the second starts after the first ends) are modeled with the first cross-link ending on the time point after this time step, in the second cross-link starting on the time point before this time step"""
            if tp_indx == 0:
                raise RuntimeError('Should not be called with tp_indx = 0')

            rr_pre = rr_dancecards[sat_indx][tp_indx_pre_t]
            rr_post = rr_dancecards[sat_indx][tp_indx_pre_t+1]

            if rr_post and t >= rr_post.release_time:
                return rr_post
            else:
                assert(t >= rr_pre.release_time)
                return rr_pre

        ########################################
        #  run stage one algorithm
        ########################################

        final_route_records = []

        visited_act_set = set()

        #  the seconds value of the very first time point
        tp_dt = time_getter_dc.get_tp_from_tp_indx(0,out_units='datetime')

        #  the main loop for the dynamic programming algorithm
        # note: have to get the generator again
        for tp_indx in time_getter_dc.get_tp_indcs ():
            if tp_indx % 10 == 0:
                print(('tp_indx: %d/%d')%(tp_indx,time_getter_dc.num_timepoints-1))

            #  nothing happens at the first time point index, because that's right at the end of the observation window
            if tp_indx == 0:
                continue


            # get time point values
            tp_last_dt = tp_dt
            tp_dt = time_getter_dc.get_tp_from_tp_indx(tp_indx,out_units='datetime')


            for sat_indx in range (self.num_sats):
                #  if we reach the end of time for this sat index, go to next satellite (i.e. were only looking for down links from observing sat at this point)
                if tp_dt > end_dt_by_sat_indx[sat_indx]:
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
                        if time_within(tp_last_dt,tp_dt,act.end):
                            xlnk_options.append(act)
                            # visited_act_set.add(act)

                ################
                #  determine for every cross-link option how much data volume the cross-link can deliver to  satellite sat_indx.

                # if sat_indx == 23:
                #     import ipdb
                #     ipdb.set_trace()


                #  each xlnk_candidate is a single xlnk window. For each of these candidates, there's a set of de-conflicted routes that move up to xlnk.data_vol amount of data volume through that xlnk window
                xlnk_candidates = []
                best_dv_seen = 0
                for xlnk in xlnk_options:
                    #  remove redundant calculations: if we've already found a cross-link option that can transport more data volume than  the maximum possible for the current cross-link, skip the current cross-link
                    if xlnk.data_vol <= best_dv_seen:
                        continue

                    # if sat_indx != 0 and (xlnk.sat_indx == 0 or xlnk.xsat_indx==0) :
                    #     import ipdb
                    #     ipdb.set_trace()
                    #  figure out which route record corresponded to the time right before this cross-link started. we can't assume it's only one time point in the past, because the cross-link could stretch across multiple timesteps
                    tp_indx_pre_xlnk = time_getter_dc.get_tp_indx_pre_t(xlnk.start,in_units='datetime')

                    # rr_last_xlnk = rr_dancecards[sat_indx][tp_indx_pre_xlnk]
                    rr_last_xlnk = get_best_rr(xlnk.start,tp_indx_pre_xlnk,sat_indx)
                    
                    #  get the route record for the corresponding crosslink partner satellite
                    xsat_indx=xlnk.get_xlnk_partner(sat_indx)
                    # rr_xsat = rr_dancecards[xsat_indx][tp_indx_pre_xlnk]
                    rr_xsat = get_best_rr(xlnk.start,tp_indx_pre_xlnk,xsat_indx)

                    #  remove redundant calculations: if the other satellite has less data volume at this time point than we do, then we can't get any more data volume from them, and it's pointless to consider it
                    if rr_last_xlnk.dv > rr_xsat.dv:
                        continue

                    #  need to figure out what data on the other satellite is data that we have not yet received on sat_indx. this returns a set of de-conflicted routes that all send valid, non-duplicated data to sat_indx
                    deconf_rts = rr_last_xlnk.get_deconflicted_routes(rr_xsat,min_dv=self.min_s_route_dv,verbose= False)

                    #todo: hmm, do we also want to consider other direction, where rr_sat is incumbent and rr_last_xlnk are the routes to deconflict?

                    #  now we need to figure out how many of these routes we can use, based upon the available crosslink bandwidth
                    xlnk_dv = xlnk.data_vol
                    xlnk_candidate_rts = []
                    # Cumulative data volume used on the crosslink by the route candidates. must be less than the available crosslink data volume
                    x_cum_dv = 0
                    # deconf_rt is a DeconflictedRoute namedtuple, from above
                    for deconf_rt in deconf_rts:
                        # if we're less than the total crosslink data volume
                        if x_cum_dv + deconf_rt.available_dv <= xlnk_dv:
                            xlnk_candidate_rts.append(deconf_rt)
                            x_cum_dv += deconf_rt.available_dv
                        # also handle the case where we don't have enough available data volume to fulfill all of the potential for route, but we can still grab enough data volume for that route ( greater than specified minimum) to use it
                        else:
                            available_dv = xlnk_dv - x_cum_dv 
                            #  check if greater than specified minimum
                            if available_dv >= self.min_s_route_dv:
                                deconf_rt.available_dv = available_dv
                                xlnk_candidate_rts.append(deconf_rt)
                                x_cum_dv += available_dv


                    #  make a new candidate entry with a record of the new data volume that we'll have if we choose that candidate.  also bring along some other relevant objects for bookkeeping
                    new_dv = rr_last_xlnk.dv + x_cum_dv
                    best_dv_seen = max(best_dv_seen,x_cum_dv)
                    # XLNK_CANDIDATES:       0         1            2                3
                    xlnk_candidates.append((new_dv,rr_last_xlnk,xlnk_candidate_rts,xlnk))

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
                    rr_new.release_time = xlnk_wind.end

                    # if sat_indx == 16:
                    #     import ipdb
                    #     ipdb.set_trace()
                    
                    # deconf_rt is a DeconflictedRoute namedtuple, from above
                    for deconf_rt in xlnk_candidate_rts:
                        #  again, make a copy because the existing data routes need to be left as is for future consideration and the algorithm
                        new_dr = copy(deconf_rt.dr)
                        #  the available data volume is simply that of the de-conflicted route, because we already made sure above that the sum of the available data volume from all de-conflicted routes is less than or equal to the throughput of the cross-link, and none of that data volume conflicts  with data volume already present in rr_new
                        new_dr.data_vol = deconf_rt.available_dv
                        new_dr.ID = dr_uid
                        #  we did not append the cross-link window to the route before ( for sake of efficiency), so now we need to do it. The second argument below records the fact that the cross-link started on the cross-link partner,  not on satellite sat_indx
                        new_dr.append_wind_to_route(xlnk_wind,window_start_sat_indx=xlnk_wind.get_xlnk_partner(sat_indx))
                        
                        rr_new.routes.append(new_dr)
                        dr_uid +=1
                else:
                    #  need to make a copy here because the last route record may be visited again, and needs to be left as is
                    rr_new = copy(rr_last)
                
                #  update the dance card at current time point.
                # IMPORTANT NOTE: this makes a new route record at the time point after the time step in which a cross-link occurred. This is the reason for using the get_best_rr() function, because an activity may start after the release time of rr_new during the same time step, and we want to consider such route records.
                rr_dancecards[sat_indx][tp_indx] = rr_new


                ################
                # go through all the activities active for this time step. pinch off any routes that can end at downlinks.
                #  wanted to wait until after performing all cross-links to check this, because cross-links may have delivered some data volume before the downlink starts in the midst of this timestep
                for act in acts:
                    # if we have already considered this activity, keep going
                    # Technically this shouldn't be necessary, because a down/crosslink should only have a start/end time within one time step. but performing this type of check here should be more efficient than checking the start time again
                    if act in visited_act_set: continue

                    if type(act) == DlnkWindow:

                        if time_within(tp_last_dt,tp_dt,act.start):

                            #  figure out if we want to use the route record from the last timepoint, or if a cross-linked has delivered more data volume. grab the latter if valid
                            rr_pre_dlnk = get_best_rr(act.start,tp_indx-1,sat_indx)

                            #  if the route record actually shows some data volume arriving at this satellite, AND we have a downlink, then we have found an optimal route to the downlink. Create a final route record and save for later.
                            if rr_pre_dlnk.dv > 0:
                                # rr_dlnk = copy(rr_pre_dlnk)
                                #  note: release time no longer matters because were not putting this route record back into the dance card.
                                rr_dlnk = RouteRecord(dv=rr_pre_dlnk.dv,release_time=act.end,routes=[])
                                
                                # if sat_indx == 16:
                                #     import ipdb
                                #     ipdb.set_trace()
                                
                                # available data volume is limited by how much we had at the last route record (sum of all data routes ending at that route record) and the downlink data volume
                                available_dv = min(rr_pre_dlnk.dv,act.data_vol)
                                rr_dlnk.dv = available_dv

                                #  now we need to update the routes within the route record to include the downlink window
                                for dr in rr_pre_dlnk.routes:
                                    dv_slice = min(dr.data_vol,available_dv)
                                    # we copied the data route objects with the copy call above, so this is not dangerous
                                    # if dv_slice < 1:
                                    #     import ipdb
                                    #     ipdb.set_trace()
                                    new_dr = copy(dr)
                                    new_dr.data_vol = dv_slice
                                    new_dr.ID = dr_uid
                                    new_dr.append_wind_to_route(act,window_start_sat_indx=sat_indx)
                                    rr_dlnk.routes.append(new_dr)

                                    dr_uid +=1
                                    available_dv -= dv_slice

                                    #  we run out of available data volume so no more routes can be grabbed
                                    if available_dv < self.dv_epsilon:
                                        break

                                # available data volume should be 0
                                assert(available_dv < self.dv_epsilon)

                                #  sanity check: also make sure that data volume for the route record is equal to the sum of the sum of the data volumes from all of the routes within it
                                assert(abs(rr_dlnk.dv - sum(dr.data_vol for dr in rr_dlnk.routes)) < self.dv_epsilon)

                                final_route_records.append(rr_dlnk)
                                visited_act_set.add(act)


        self.final_route_records = final_route_records

        all_routes = [dr for rr in final_route_records for dr in rr.routes]
        for dr in all_routes:
            dr.validate_route()

        return all_routes, dr_uid


    def get_stats(self,verbose=False):
        stats = {}
        # stats['num_dlnk_windows'] = sum([len (self.dlnk_winds_flat[sat_indx]) for sat_indx in range (self.num_sats)])
        # stats['num_xlnk_windows'] = sum([len (self.xlnk_winds[sat_indx][xsat_indx]) for xsat_indx in range (self.num_sats) for sat_indx in range (self.num_sats)])

        if verbose:
            for rr_indx, rr in  enumerate (self.final_route_records):
                print ('route record %d, %d routes, sum rts/ dlnk dv: %.0f/%.0f Mb (%.2f%%)'%(rr_indx,len(rr.routes),sum(dr.data_vol for dr in rr.routes),rr.routes[0].get_dlnk().data_vol,100*sum(dr.data_vol for dr in rr.routes)/rr.routes[0].get_dlnk().data_vol))
                for dr_indx, dr in  enumerate (rr.routes):
                    pass
                    # print('  dr %d %s'%(dr_indx,dr))
                    print('  dr %d %s'%(dr_indx,dr.get_dlnk ()))

            # print ( "Obs dv: %f" % ( self.obs_wind.data_vol))
            # print ( "Considering %d downlink windows" % (stats['num_dlnk_windows']))
            # print ( "Considering %d crosslink windows" % (stats['num_xlnk_windows']))

        return stats




