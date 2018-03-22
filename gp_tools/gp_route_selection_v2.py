# Algorithm for creating a set of data routes for a given observation, later fed into the activity scheduling stage
# 
# This version uses a dynamic programming algorithm to optimize selected paths
# 
# @author Kit Kennedy
#
#  note that a path is the same as a route. 

from  datetime import timedelta
from copy import  deepcopy

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

class GPDataRouteSelection():
    """docstring for GP route selection"""

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
        rs_general_params = gp_params['gp_general_params']['route_selection_general_params']
        rs_params = gp_params['gp_general_params']['route_selection_params_v2']
        gp_general_other_params = gp_params['gp_general_params']['other_params']
        gp_orbit_link_params = gp_params['gp_orbit_link_params']
        gp_inst_params = gp_params['gp_instance_params']['route_selection_params']

        self.num_sats=sat_params['num_sats']
        self.num_paths=rs_params['num_paths']
        #  these times indicate the (largest) window over which we are considering routes
        self.sel_start_utc_dt  = tt.iso_string_to_dt (gp_inst_params['start_utc'])
        self.sel_end_utc_dt  = tt.iso_string_to_dt (gp_inst_params['end_utc'])

        # get the smallest time step used in orbit link. this is the smallest time step we need to worry about in data route selection
        self.act_timestep = min(gp_orbit_link_params['xlnk_max_len_s'],gp_orbit_link_params['dlnk_max_len_s'])

        self.min_path_dv =rs_params['min_path_dv_Mb']
        self.wind_filter_duration =  timedelta (seconds =rs_general_params['wind_filter_duration_s'])
        self.latency_params =  gp_general_other_params['latency_calculation']

        if self.latency_params['obs'] not in ['start','end']:
            raise NotImplementedError
        if self.latency_params['dlnk'] not in ['start','end','center']:
            raise NotImplementedError


    @staticmethod
    def  filter_windows(obs_wind,dlink_winds_flat,xlink_winds,num_sats,end_utc_dt,wind_filter_duration):
        first =  obs_wind.end
        last =  min ( end_utc_dt, first +  wind_filter_duration)

        dlink_winds_flat_filtered = [[] for sat_indx in  range (num_sats)]
        xlink_winds_flat_filtered = [[[] for xsat_indx in  range ( num_sats)] for sat_indx in  range (num_sats)]

        for sat_indx in  range (num_sats):
            for xsat_indx in  range ( num_sats):
                for wind in xlink_winds[sat_indx][xsat_indx]:
                    if  wind.start > first  and  wind.end  <last:
                        xlink_winds_flat_filtered[sat_indx][xsat_indx]. append ( wind)

            for wind in dlink_winds_flat[sat_indx]:
                if  wind.start > first  and  wind.end  <last:
                    dlink_winds_flat_filtered[sat_indx]. append ( wind)

        return dlink_winds_flat_filtered, xlink_winds_flat_filtered

    def run ( self,obs_wind,dlink_winds_flat,xlink_winds, verbose = True):

        # construct a set of dance cards for every satellite, 
        # each of which keeps track of all of the activities of satellite 
        # can possibly execute at any given time slice delta T. 
        #  start the dance card from the end of the observation window, because we are only selecting routes after the observation
        act_dancecards = [Dancecard(obs_wind.end,self.sel_end_utc_dt,self.act_timestep) for sat_indx in range (self.num_sats)]
        
        #  make another set of dance cards that keep track of the path data structure at a given time point
        path_dancecards = [Dancecard(obs_wind.end,self.sel_end_utc_dt,self.act_timestep,item_init=None) for sat_indx in range (self.num_sats)]

        for sat_indx in range (self.num_sats): 
            act_dancecards[sat_indx].add_winds_to_dancecard(dlink_winds_flat[sat_indx])
            act_dancecards[sat_indx].add_winds_to_dancecard(xlink_winds[sat_indx])





    def get_stats(self,verbose=True):
        stats = {}
        stats['num_dlnk_windows'] = sum([len (self.dlink_winds_flat[sat_indx]) for sat_indx in range (self.num_sats)])
        stats['num_xlnk_windows'] = sum([len (self.xlink_winds[sat_indx][xsat_indx]) for xsat_indx in range (self.num_sats) for sat_indx in range (self.num_sats)])

        if verbose:
            print ( "Obs dv: %f" % ( self.obs_wind.data_vol))
            print ( "Considering %d downlink windows" % (stats['num_dlnk_windows']))
            print ( "Considering %d crosslink windows" % (stats['num_xlnk_windows']))

        return stats




