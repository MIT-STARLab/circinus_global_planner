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

    def get_deconflicted_routes(self,rr_other,min_dv=0):

        deconflicted_routes = []

        #  if self has routes, then it's possible that these routes share some of the same data from the routes in the other route record ( on the other satellite). this could come up because the routes share some of the same initial windows. we need to determine where the routes split, and use that to determine how much data volume is actually unique ( hasn't already arrived on self), and available to receive
        if len(self.routes) > 0:
            for dr_1 in self.routes:
                for dr_2 in rr_other.routes:
                    #  get the last window the routes had in common
                    split_wind = dr_1.get_split(dr_2)
                    available_dv = split_wind.data_vol - dr_1.data_vol

                    if available_dv >= min_dv and no_route_conflict(dr_1,dr_2): 
                        deconflicted_routes.append(DeconflictedRoute(dr=dr_2,available_dv=available_dv))

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

        self.min_path_dv =rs_general_params['min_path_dv_Mb']
        #  the "effectively zero" number. defining in terms of the minimum path data volume for future proofing
        self.dv_epsilon = self.min_path_dv / 1000
        self.wind_filter_duration =  timedelta (seconds =rs_general_params['wind_filter_duration_s'])
        self.latency_params =  gp_general_other_params['latency_calculation']

        if self.latency_params['obs'] not in ['start','end']:
            raise NotImplementedError
        if self.latency_params['dlnk'] not in ['start','end','center']:
            raise NotImplementedError


    @staticmethod
    def  filter_windows(obs_wind,dlnk_winds_flat,xlnk_winds,num_sats,end_utc_dt,wind_filter_duration):
        first =  obs_wind.end
        last =  min ( end_utc_dt, first +  wind_filter_duration)

        dlink_winds_flat_filtered = [[] for sat_indx in  range (num_sats)]
        xlink_winds_flat_filtered = [[[] for xsat_indx in  range ( num_sats)] for sat_indx in  range (num_sats)]

        for sat_indx in  range (num_sats):
            for xsat_indx in  range ( num_sats):
                for wind in xlnk_winds[sat_indx][xsat_indx]:
                    if  wind.start > first  and  wind.end  <last:
                        xlink_winds_flat_filtered[sat_indx][xsat_indx]. append ( wind)

            for wind in dlnk_winds_flat[sat_indx]:
                if  wind.start > first  and  wind.end  <last:
                    dlink_winds_flat_filtered[sat_indx]. append ( wind)

        return dlink_winds_flat_filtered, xlink_winds_flat_filtered

    def run ( self,obs_wind,dlnk_winds_flat,xlnk_winds, dr_uid, verbose = False):
        #  note that a potential source of slowness in the code below is the creation of new RouteRecord objects for every sat at every time step 
        # TODO: figure out if there's a more efficient way to do this -  for example, we shouldn't have to preserve these objects when they're more than a certain number of time steps in the past


        start_dt = obs_wind.end
        end_dt = min(self.sel_latest_end_dt,start_dt + self.wind_filter_duration)

        dlnk_winds_flat_filt,xlnk_winds_filt =  self.filter_windows (obs_wind,dlnk_winds_flat,xlnk_winds, self.num_sats, self.sel_latest_end_dt, self.wind_filter_duration)

        print ('Running route selection for obs: %s'%(obs_wind))
        print ('from %s to %s'%(start_dt, end_dt))

        # construct a set of dance cards for every satellite, 
        # each of which keeps track of all of the activities of satellite 
        # can possibly execute at any given time slice delta T. 
        #  start the dance card from the end of the observation window, because we are only selecting routes after the observation
        act_dancecards = [Dancecard(start_dt,end_dt,self.act_timestep,mode='timestep') for sat_indx in range (self.num_sats)]
        
        #  make another set of dance cards that keep track of the route record data structures at a given time point
        #  note that even though these dance cards are in different modes, they share the same time points and time steps (both values and indices)
        rr_dancecards = [Dancecard(start_dt,end_dt,self.act_timestep,item_init=None,mode='timepoint') for sat_indx in range (self.num_sats)]

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

        #  all dancecards here share the same time point indices
        time_getter_dc = act_dancecards[0]
        timepoint_indcs = time_getter_dc.get_tp_indcs ()

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
        for tp_indx in timepoint_indcs:
            if tp_indx % 10 == 0:
                print(('tp_indx: %d/%d')%(tp_indx,len(timepoint_indcs)-1))

            #  nothing happens at the first time point index, because that's right at the end of the observation window
            if tp_indx == 0:
                continue


            # get time point values
            tp_last_dt = tp_dt
            tp_dt = time_getter_dc.get_tp_from_tp_indx(tp_indx,out_units='datetime')

            for sat_indx in range (self.num_sats):
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
                    deconf_rts = rr_last_xlnk.get_deconflicted_routes(rr_xsat,min_dv=self.min_path_dv)

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
                            if available_dv >= self.min_path_dv:
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
                        #  we did not append the cross-link window to the route before ( for sake of efficiency), so now we need to do it. The second argument below records the fact that the cross-link started on the cross-link partner,  not on satellite sat_indx
                        new_dr.data_vol = deconf_rt.available_dv
                        new_dr.ID = dr_uid
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
        return all_routes


    def get_stats(self,verbose=False):
        stats = {}
        # stats['num_dlnk_windows'] = sum([len (self.dlnk_winds_flat[sat_indx]) for sat_indx in range (self.num_sats)])
        # stats['num_xlnk_windows'] = sum([len (self.xlnk_winds[sat_indx][xsat_indx]) for xsat_indx in range (self.num_sats) for sat_indx in range (self.num_sats)])

        if verbose:
            for rr_indx, rr in  enumerate (self.final_route_records):
                print ('route record %d, %d routes, sum rts/ dlnk dv: %f/%f Mb '%(rr_indx,len(rr.routes),sum(dr.data_vol for dr in rr.routes),rr.routes[0].get_dlnk().data_vol))
                for dr_indx, dr in  enumerate (rr.routes):
                    print('  %d %s'%(dr_indx,dr))

            # print ( "Obs dv: %f" % ( self.obs_wind.data_vol))
            # print ( "Considering %d downlink windows" % (stats['num_dlnk_windows']))
            # print ( "Considering %d crosslink windows" % (stats['num_xlnk_windows']))

        return stats




