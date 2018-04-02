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

    # multi:  there are forks and convergences in the route; data can split out from one window to travel through multiple windows and converge back on another window.  temporal consistency must still be maintained though

class DataMultiRoute():

    def __init__(self,ID,data_routes):
        self.ID = ID
        self.data_routes = data_routes
        self.data_vol_by_dr = {dr:dr.data_vol for dr in data_routes}

        for dr in data_routes:
            if not type(dr) == DataRoute:
                raise RuntimeError('only data route objects should be used to construct a new data multi-route')

    @property
    def data_vol(self):
        return sum(self.data_vol_by_dr.values())

    def get_winds(self):
        return (wind for dr in self.data_routes for wind in dr.get_winds())

    def get_obs( self):
        #  use first route because all routes should have the same observation
        return self.data_routes[0].get_obs()

    def get_dlnk(self):
        #  use first route because all routes should have the same downlink
        return self.data_routes[0].get_dlnk()

    def has_xlnk(self):
        return any(dr.has_xlnk() for dr in self.data_routes)

    def __repr__(self):
        return  '(DataMultiRoute %d, routes: %s)'%(self.ID,{dr.ID: self.data_vol_by_dr[dr] for dr in self.data_routes})

    def accumulate_dr( self, candidate_dr,min_dmr_candidate_dv=0):
        """[summary]
        
        [description]
        :param dmrs: data multi routes to combine into a new multiroute
        :type dmrs: list(DataMultiRoute)
        """

        # need to have matching observation and downlink for candidate to be added on to the multi-route
        if not candidate_dr.get_obs() == self.get_obs():
            return False
        if not candidate_dr.get_dlnk() == self.get_dlnk():
            return False

        avail_dv_by_wind = {}
        # req_dv_by_wind = {}
        for dr in self.data_routes:
            for wind in dr.route:
                #  if we didn't yet encounter this window in any of the routes in self
                avail_dv_by_wind.setdefault(wind,wind.data_vol)
                avail_dv_by_wind[wind] -= dr.data_vol
                # req_dv_by_wind.setdefault(wind,0)
                # req_dv_by_wind[wind] += dr.data_vol


                # if avail_check >= dr.data_vol:
                #     avail_dv_by_wind[wind] -= dr.data_vol
                #     # record the fact that this window is contained in dr_2
                #     other_rts_by_wind.setdefault(wind,[]).append(dr_2_indx)

        # candidate_feasible = True
        candidate_dv =candidate_dr.data_vol
        for wind in candidate_dr.route:
            # if wind in avail_dv_by_wind.keys():
            #     new_dv = avail_dv_by_wind[wind] - candidate_dr.data_vol
            usable_wind_dv = min(avail_dv_by_wind.get(wind,wind.data_vol),candidate_dr.data_vol)

            candidate_dv = min(candidate_dv, usable_wind_dv)

        if candidate_dv > min_dmr_candidate_dv:
            self.data_routes.append(candidate_dr)
            self.data_vol_by_dr[candidate_dr] = candidate_dv
            return True
        else:
            return False


class DataRoute():
    '''
    Contains all the relevant information about the path taken by a single data packet traveling through the satellite network
    '''

    # note this route is simple:  there are no forks in the route; there is a simple linear path from an observation to a downlink through which data flows. all windows must be in temporal order.

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

    def get_winds(self):
        return (wind for wind in self.route)

    def get_obs( self):
        return self.route[0]

    def get_dlnk( self):
        return self.route[-1]

    def has_xlnk(self):
        for wind in self.route:
            if type(wind) == XlnkWindow:
                return True

        return False

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

                #  if this happens to be a unidirectional window, and the current satellite index is not the transmitting satellite for that window, there's a problem
                if not wind.symmetric and curr_sat_indx != wind.tx_sat:
                    string ='routing_objects.py: Found incorrect tx sat at window indx %d in route. Route string: %s'%( windex, self.get_route_string())
                    raise RuntimeError(string)

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

    def count_overlap(self,other,window_option ='shared_window'):

        overlap_count = 0
        for wind1 in self.route:
            if type(wind1) != XlnkWindow:
                continue

            for wind2 in  other.route:
                if type(wind2) != XlnkWindow:
                    continue

                if wind1 == wind2:
                    if window_option =='shared_window':
                        #  only count true overlaps, meaning there's not enough space in the window for both routes
                        if self.data_vol + other.data_vol > wind1.data_vol: 
                            overlap_count += 1
                    elif window_option == 'mutex_window':
                        overlap_count += 1
                    else:
                        raise NotImplementedError


        return overlap_count

    def is_overlapping(self,other,window_option ='shared_window'):

        overlap_count = 0
        for wind1 in self.route:
            if type(wind1) != XlnkWindow:
                continue

            for wind2 in  other.route:
                if type(wind2) != XlnkWindow:
                    continue

                if wind1 == wind2:
                    if window_option =='shared_window':
                        #  only count true overlaps, meaning there's not enough space in the window for both routes
                        if self.data_vol + other.data_vol > wind1.data_vol: 
                            return True
                    elif window_option == 'mutex_window':
                        return True
                    else:
                        raise NotImplementedError

        return False



class LinkInfo(object):
    """docstring fos  LinkInfo"""
    def __init__(self,data_routes=[],total_data_vol=0,used_data_vol=0):
        self.data_routes = data_routes
        self.total_data_vol = total_data_vol
        self.used_data_vol = used_data_vol

    def __str__( self):
        return  "routes: "+str(self.data_routes) + " ; dv %.0f/%.0f Mb" % ( self.used_data_vol, self.total_data_vol)


