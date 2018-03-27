# Data structures for use in routing planning code
# 
# @author Kit Kennedy

from copy import copy

from circinus_tools  import  constants as const
from .custom_activity_window import   ObsWindow,  DlnkWindow, XlnkWindow

DATE_STRING_FORMAT = 'short'
# DATE_STRING_FORMAT = 'iso'

def short_date_string(dt):
    return dt.strftime("%H:%M:%S")

def date_string(dt):
    if DATE_STRING_FORMAT == 'iso':
        return dt.isoformat()
    if DATE_STRING_FORMAT == 'short':
        return  short_date_string(dt)

class DataRoute(object):
    '''
    Contains all the relevant information about the path taken by a single data packet traveling through the satellite network
    '''

    def __init__(self, ID, route  =[], window_start_sats={},dv=0):

        self.ID =  ID

        # the list storing all objects in the route; a list ObsWindow, XlnkWindow, XlnkWindow...DlnkWindow
        self.route =  route

        # this keeps track, for each window along the route, of the satellite (sat_indx) that the data was on at the beginning of the window. This is necessary because the window objects themselves are used across paths and do not store information about which sense they are used in
        #  dictionary with keys being the window objects themselves
        self.window_start_sats = window_start_sats

        self.data_vol = dv

        #  the amount of capacity on the path that actually ends up scheduled for usage
        self.scheduled_dv = const.UNASSIGNED

        self.sort_windows()   

        #  check the timing and satellite indices along the route.  throws an exception if a problem is seen
        self.validate_route()

    def __copy__(self):
        newone = type(self)(self.ID,dv=self.data_vol)
        #  make a shallow copy of these container objects -  we want to refer to the same nested objects within the containers, but want a new container in both cases
        newone.route = copy(self.route)
        newone.window_start_sats = copy(self.window_start_sats)
        return newone

    def append_wind_to_route( self,wind,window_start_sat_indx):
        self.route.append(wind)
        self.window_start_sats[wind] = window_start_sat_indx

    def get_obs( self):
        return self.route[0]

    def get_dlnk( self):
        return self.route[-1]

    def get_latency( self,units='minutes',obs_option = 'end', dlnk_option = 'center'):
        obs =  self.route[0]
        dlnk =  self.route[-1]

        lat_start = getattr(obs,obs_option)
        lat_end = getattr(dlnk,dlnk_option)
    
        if units == 'minutes':
            return (lat_end-lat_start).total_seconds()/60
        else:
            raise NotImplementedError

    def  sort_windows(self):
        self.route.sort(key=lambda x: x.start)

        # #  do some basic validation that the windows are in order
        # last_wind_end = const.UNASSIGNED_DT_NEG_INF
        # for wind in self.route:
        #     if wind.start <last_wind_end  or wind.end  < last_wind_end:
        #         raise Exception ('routing_objects.py:  windows are out of order')
        #     last_wind_end = wind.end

    def  get_route_string( self,  time_base= None):
        out_string = "dr %d: "% ( self.ID)

        for wind in self.route:

            if type (wind)  == ObsWindow:
                start_str =  "%.0fs" % ( wind.start-time_base).total_seconds() if  time_base else  date_string(wind.start)
                end_str =  "%.0fs" % ( wind.end-time_base).total_seconds() if  time_base else  date_string(wind.end)
                out_string  +=  "o %d dv %.0f %s,%s" % (wind.sat_indx, wind.data_vol,start_str,end_str)
            elif type (wind)  == XlnkWindow:
                start_str =  "%.0fs" % ( wind.start-time_base).total_seconds() if  time_base else  date_string(wind.start)
                end_str =  "%.0fs" % ( wind.end-time_base).total_seconds() if  time_base else  date_string(wind.end)
                sat_indx=self.window_start_sats[wind]
                # xsat_indx=wind.xsat_indx  if self.window_start_sats[wind] == wind.sat_indx else wind.sat_indx
                xsat_indx=wind.get_xlnk_partner(self.window_start_sats[wind])
                out_string  +=  " -> x %d,%d dv %.0f %s,%s" % (sat_indx, xsat_indx, wind.data_vol,start_str,end_str)
            elif type (wind)  == DlnkWindow:
                start_str =  "%.0fs" % ( wind.start-time_base).total_seconds() if  time_base else  date_string(wind.start)
                end_str =  "%.0fs" % ( wind.end-time_base).total_seconds() if  time_base else  date_string(wind.end)
                sat_indx= wind.sat_indx
                out_string  +=  " -> d %d dv %.0f %s,%s" % (sat_indx, wind.data_vol,start_str,end_str)
        
        return out_string

    def __repr__(self):
        if not self.scheduled_dv:
            return  '('+self.get_route_string()+')'
        else:
            if self.scheduled_dv == const.UNASSIGNED:
                return  '('+self.get_route_string()+'; sched dv: %s/%.0f Mb'%( 'none', self.data_vol)+')'
            else:
                return  '('+self.get_route_string()+'; sched dv: %.0f/%.0f Mb'%( self.scheduled_dv, self.data_vol)+')'

    def __getitem__(self, key):
        """ getter for internal route by index"""
        return self.route[key]

    def validate_route (self):
        """ validates timing and ordering of route
        
        [description]
        :raises: Exception, Exception, Exception
        """

        if len( self.route) == 0:
            return

        obs = self.route[0]
        if not type (obs) is ObsWindow:
            raise Exception('First window on route was not an ObsWindow instance. Route string: %s'%( self.get_route_string()))

        curr_sat_indx = obs.sat_indx
        next_sat_indx = obs.sat_indx
        last_time = obs.start

        #  trace through the route and make sure: 1. we cross through satellites in order and 2.  every activity along the path starts after the last activity ended
        for windex, wind in  enumerate(self.route):
            if self.window_start_sats[wind] != next_sat_indx:
                string = 'routing_objects.py: Found the incorrect sat indx at window indx %d in route. Route string: %s'%( windex, self.get_route_string())
                raise RuntimeError(string)
            if not wind.start >= last_time or not wind.end >= last_time:
                string ='routing_objects.py: Found a bad start time at window indx %d in route. Route string: %s'%( windex, self.get_route_string())
                raise RuntimeError( string)

            if not self.data_vol <= wind.data_vol:
                string ='routing_objects.py: Found bad dv at window indx %d in route. Route string: %s'%( windex, self.get_route_string())
                raise RuntimeError( string)

            #  note that we manually trace the satellite index through cross-link window here. This is a bit redundant with the functionality of window_start_sats,  but adds a little bit more of a warm, happy, comfortable feeling in the area checking
            if type (wind) is XlnkWindow:
                curr_sat_indx = next_sat_indx
                next_sat_indx=wind.get_xlnk_partner(curr_sat_indx)
                # next_sat_indx = wind.xsat_indx if not wind.xsat_indx == curr_sat_indx else wind.sat_indx

            last_time = wind.end

    def get_split(self,other):
        """ return the last window in common with other data route
        
        iterate through the windows in both routes in temporal order and return the last window that is the same for the two routes
        :param other:  another data route
        :type other: DataRoute
        :returns:  last common window
        :rtype: {ActivityWindow}
        """

        #  route should always start the same place, the initial observation window
        assert self.route[0] == other.route[0]

        len_other = len(other.route)

        #  this will record the index of the last common window of the two routes
        #  increment every time a window is proven to be contained in both routes
        split_windex = -1
        for windex,wind in enumerate(self.route):

            #  if we reached a window in self that is longer than the route for other
            if windex+1 > len_other:
                break

            #  note that this tests the window ID for equality
            if wind != other.route[windex]:
                break

            #  if we've made it here, then both of the data routes must have this window in common
            split_windex += 1

        return self.route[split_windex]

    def count_overlap(self,other):


        overlap_count = 0
        for wind1 in self.route:
            if type(wind1) != XlnkWindow:
                continue

            for wind2 in  other.route:
                if type(wind2) != XlnkWindow:
                    continue

                if wind1 == wind2:
                    overlap_count += 1


        return overlap_count



class LinkInfo(object):
    """docstring fos  LinkInfo"""
    def __init__(self,data_routes=[],total_data_vol=0,used_data_vol=0):
        self.data_routes = data_routes
        self.total_data_vol = total_data_vol
        self.used_data_vol = used_data_vol

    def __str__( self):
        return  "routes: "+str(self.data_routes) + " ; dv %.0f/%.0f Mb" % ( self.used_data_vol, self.total_data_vol)


