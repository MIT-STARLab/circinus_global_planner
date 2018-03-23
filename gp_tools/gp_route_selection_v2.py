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
from .routing_objects import DataRoute


class RouteRecord():
    """docstring for RouteRecord"""
    def __init__(self,dv,routes=[]):
        """create a new route record
        
        used in dynamic programming algorithm to keep track of data volume and data routes leading to a given point in time on a satellite
        :param dv:  data volume for set of routes
        :type dv:  float
        :param routes:  list of DataRoute objects that currently end at a given time on a satellite , defaults to []
        :type routes: list, optional
        """
        self.dv = dv
        self.routes = routes

    DeconflictedRoute = namedtuple('DeconflictedRoute', 'dr available_dv')

    def get_deconflicted_routes(self,rr_other,min_dv=0):

        deconflicted_routes = []

        #  if self has routes, then it's possible that these routes share some of the same data from the routes in the other route record ( on the other satellite). this could come up because the routes share some of the same initial windows. we need to determine where the routes split, and use that to determine how much data volume is actually unique ( hasn't already arrived on self), and available to receive
        if len(self.routes) > 0:
            for dr_1 in self.routes:
                for dr_2 in rr_other.routes:
                    #  get the last window the routes had in common
                    split_wind = dr1.get_split(dr_2)
                    available_dv = split_wind.data_vol - dr_1.data_vol

                    if available_dv >= min_dv: 
                        deconflicted_routes.append(DeconflictedRoute(dr=dr_2,available_dv=available_dv))

        #  if self has no routes ( it hasn't yet received data from anywhere), then every route in  the other route record (on the other satellite) is valid
        else:
            for dr_2 in rr_other.routes:
                available_dv = dr_2.data_vol
                if available_dv >= min_dv: 
                    DeconflictedRoutected_routes.append(DeconflictedRoute(dr=dr_2,available_dv=available_dv)) 

        return deconflicted_routes

    def __copy__(self):
        newone = type(self)( self.dv)
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
        gp_inst_params = gp_params['gp_instance_params']['route_selection_params']

        self.num_sats=sat_params['num_sats']
        #  these times indicate the (largest) window over which we are considering routes
        self.sel_start_utc_dt  = tt.iso_string_to_dt (gp_inst_params['start_utc'])
        self.sel_end_utc_dt  = tt.iso_string_to_dt (gp_inst_params['end_utc'])

        # get the smallest time step used in orbit link. this is the smallest time step we need to worry about in data route selection
        self.act_timestep = min(link_params['xlnk_max_len_s'],link_params['dlnk_max_len_s'])

        self.min_path_dv =rs_general_params['min_path_dv_Mb']
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

        # construct a set of dance cards for every satellite, 
        # each of which keeps track of all of the activities of satellite 
        # can possibly execute at any given time slice delta T. 
        #  start the dance card from the end of the observation window, because we are only selecting routes after the observation
        act_dancecards = [Dancecard(obs_wind.end,self.sel_end_utc_dt,self.act_timestep,mode='timestep') for sat_indx in range (self.num_sats)]
        
        #  make another set of dance cards that keep track of the route record data structures at a given time point
        #  note that even though these dance cards are in different modes, they share the same time points and time steps (both values and indices)
        rr_dancecards = [Dancecard(obs_wind.end,self.sel_end_utc_dt,self.act_timestep,item_init=None,mode='timepoint') for sat_indx in range (self.num_sats)]

        obs_dv = obs_wind.data_vol

        for sat_indx in range (self.num_sats): 
            act_dancecards[sat_indx].add_winds_to_dancecard(dlnk_winds_flat[sat_indx])
            #  cross-link Windows matrix is symmetric ( and upper triangular)
            for xsat_indx in  range (sat_indx+1, self.num_sats):
                act_dancecards[sat_indx].add_winds_to_dancecard(xlnk_winds[sat_indx][xsat_indx])
            #  first point for every non-observing satellite has a zero data volume, no path object
            rr_dancecards[sat_indx][0] = RouteRecord(dv=0,routes=[])

        #  observing sat starts with an initial data route that includes only the observation window
        first_dr = DataRoute(dr_uid, route =[obs_wind], window_start_sats={obs_wind: obs_wind.sat_indx},dv=obs_dv)
        #  put this initial data route on the first route record for the observing satellite
        rr_dancecards[obs_wind.sat_indx][0] = RouteRecord(dv=obs_dv,routes=[first_dr])

        #  all dancecards here share the same time point indices
        time_getter_dc = act_dancecards[0]
        timepoint_indcs = time_getter_dc.get_tp_indcs ()

        def time_within(t1,t2,toi):
            return toi >= t1 and toi < t2


        ########################################
        #  run stage one algorithm
        ########################################

        final_route_records = []

        visited_act_set = set()

        #  the seconds value of the very first time point
        tp_dt = time_getter_dc.get_tp_from_tp_indx(0,out_units='datetime')

        #  the main loop for the dynamic programming algorithm
        for tp_indx in timepoint_indcs:
            #  nothing happens at the first time point index, because that's right at the end of the observation window
            if tp_indx == 0:
                continue

            # get time point values
            tp_last_dt = tp_dt
            tp_dt = time_getter_dc.get_tp_from_tp_indx(tp_indx,out_units='datetime')

            for sat_indx in range (self.num_sats):
                #  get the activities that were active during the time step immediately preceding time point
                acts = act_dancecards[sat_indx].get_objects_pre_tp_indx(tp_indx)

                #  get last route record
                rr_last = rr_dancecards[sat_indx][tp_indx-1]

                xlnk_options =[]

                ################
                # go through all the activities active for this time step. pinch off any routes that can end at downlinks. make a list of cross-link windows that are candidates for sending data to sat_indx
                for act in acts:
                    # if we have already considered this activity, keep going
                    # Technically this shouldn't be necessary, because a down/crosslink should only have a start/end time within one time step. but performing this type of check here should be more efficient than checking the start time again
                    if act in visited_act_set: continue

                    if type(act) == DlnkWindow:

                        #  if the last route record actually shows some data volume arriving at this satellite, AND we have a downlink, then we have found an optimal route to the downlink. Create a final route record and save for later.
                        if rr_last.dv > 0 and time_within(tp_last_dt,tp_dt,act.start):
                            # import ipdb
                            # ipdb.set_trace()
                            rr_dlnk = copy(rr_last)
                            # available data volume is limited by how much we had at the last route record (sum of all data routes ending at that route record) and the downlink data volume
                            available_dv = min(rr_last.dv,act.data_vol)
                            rr_dlnk.dv = available_dv

                            #  now we need to update the routes within the route record to include the downlink window
                            for dr in rr_dlnk.routes:
                                dv_slice = min(dr.data_vol,available_dv)
                                # we copied the data route objects with the copy call above, so this is not dangerous
                                dr.data_vol = dv_slice
                                dr.ID = dr_uid
                                dr.append_wind_to_route(act,window_start_sat_indx=sat_indx)

                                dr_uid +=1
                                available_dv -= dv_slice

                            final_route_records.append(rr_dlnk)
                            visited_act_set.add(act)

                    if type(act) == XlnkWindow:
                        if time_within(tp_last_dt,tp_dt,act.end):
                            xlnk_options.append(act)
                            visited_act_set.add(act)

                ################
                #  determine for every cross-link option how much data volume the cross-link can deliver to  satellite sat_indx.
                xlnk_candidates = []
                for xlnk in xlnk_options:
                    #  figure out which route record corresponded to the time right before this cross-link started. we can't assume it's only one time point in the past, because the cross-link could stretch across multiple timesteps
                    tp_indx_pre_xlnk = time_getter_dc.get_tp_indx_pre_t(xlnk.start,in_units='datetime')
                    rr_last_xlnk = rr_dancecards[sat_indx][tp_indx_pre_xlnk]
                    
                    #  get the route record for the corresponding crosslink partner satellite
                    xsat_indx=xlnk.get_xlnk_partner(sat_indx)
                    rr_xsat = rr_dancecards[xsat_indx][tp_indx_pre_xlnk]

                    #  need to figure out what data on the other satellite is data that we have not yet received on sat_indx. this returns a set of de-conflicted routes that all send valid, non-duplicated data to sat_indx
                    deconf_rts = rr_last_xlnk.get_deconflicted_routes(rr_xsat,min_dv=self.min_path_dv)

                    #  now we need to figure out how many of these routes we can use, based upon the available crosslink bandwidth
                    xlnk_dv = xlnk.data_vol
                    xlnk_candidate_rts = []
                    x_cum_dv = 0
                    # deconf_rt is a DeconflictedRoute namedtuple, from above
                    for deconf_rt in deconf_rts:
                        if x_cum_dv + deconf_rt.available_dv <= xlnk_dv:
                            xlnk_candidate_rts.append(deconf_rt)
                            x_cum_dv += deconf_rt.available_dv

                    #  make a new candidate entry with a record of the new data volume that we'll have if we choose that candidate.  also bring along some other relevant objects for bookkeeping
                    new_dv = rr_last_xlnk.dv + x_cum_dv
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
                    #  need to make a copy here because the last route record may be visited again, and needs to be left as is
                    rr_new = copy(best_xlnk_cand[1])
                    rr_new.dv = best_xlnk_cand[0]
                    xlnk_candidate_rts = best_xlnk_cand[2]
                    xlnk_wind = best_xlnk_cand[3]
                    
                    # deconf_rt is a DeconflictedRoute namedtuple, from above
                    for deconf_rt in xlnk_candidate_rts:
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
                
                rr_dancecards[sat_indx][tp_indx] = rr_new


        self.final_route_records = final_route_records

        all_routes = [dr for rr in final_route_records for dr in rr.routes]
        return all_routes


    def get_stats(self,verbose=False):
        stats = {}
        # stats['num_dlnk_windows'] = sum([len (self.dlnk_winds_flat[sat_indx]) for sat_indx in range (self.num_sats)])
        # stats['num_xlnk_windows'] = sum([len (self.xlnk_winds[sat_indx][xsat_indx]) for xsat_indx in range (self.num_sats) for sat_indx in range (self.num_sats)])

        if verbose:
            for rr_indx, rr in  enumerate (self.final_route_records):
                print ('route record %d, %d routes'%(rr_indx,len(rr.routes)))
                for dr_indx, dr in  enumerate (rr.routes):
                    print('  %d %s'%(dr_indx,dr))

            # print ( "Obs dv: %f" % ( self.obs_wind.data_vol))
            # print ( "Considering %d downlink windows" % (stats['num_dlnk_windows']))
            # print ( "Considering %d crosslink windows" % (stats['num_xlnk_windows']))

        return stats




