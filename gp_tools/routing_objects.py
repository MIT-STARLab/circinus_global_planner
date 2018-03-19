from circinus_tools  import  constants as const
from .custom_activity_window import   ObsWindow,  DlnkWindow, XlnkWindow

class DataRoute(object):
    '''
    Contains all the relevant information about the path taken by a single data packet traveling through the satellite network
    '''

    def __init__(self, ID, route  =[], window_start_sats={},dv=0):

        self.ID =  ID

        # the list storing all objects in the route; a list ObsWindow, XlnkWindow, XlnkWindow...DlnkWindow
        self.route =  route

        # this keeps track, for each window along the route, of the satellite that the data was on at the beginning of the window. This is necessary because the window objects themselves are used across paths and do not store information about which sense they are used in
        #  dictionary with keys being the window objects themselves
        self.window_start_sats = window_start_sats

        self.data_vol = dv

        #  the amount of capacity on the path that actually ends up scheduled for usage
        self.scheduled_dv = const.UNASSIGNED

        self.sort_windows()   

        #  check the timing and satellite indices along the route.  throws an exception if a problem is seen
        self.validate_route()

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
                start_str =  "%.0fs" % ( wind.start-time_base).total_seconds() if  time_base else  wind.start.isoformat()
                end_str =  "%.0fs" % ( wind.end-time_base).total_seconds() if  time_base else  wind.end.isoformat()
                out_string  +=  "o %d; %s,%s" % (wind.sat_indx,start_str,end_str)
            elif type (wind)  == XlnkWindow:
                start_str =  "%.0fs" % ( wind.start-time_base).total_seconds() if  time_base else  wind.start.isoformat()
                end_str =  "%.0fs" % ( wind.end-time_base).total_seconds() if  time_base else  wind.end.isoformat()
                sat_indx=self.window_start_sats[wind]
                xsat_indx=wind.xsat_indx  if self.window_start_sats[wind] == wind.sat_indx else wind.sat_indx
                out_string  +=  " -> x %d,%d; %s,%s" % (sat_indx, xsat_indx,start_str,end_str)
            elif type (wind)  == DlnkWindow:
                start_str =  "%.0fs" % ( wind.start-time_base).total_seconds() if  time_base else  wind.start.isoformat()
                end_str =  "%.0fs" % ( wind.end-time_base).total_seconds() if  time_base else  wind.end.isoformat()
                sat_indx= wind.sat_indx
                out_string  +=  " -> d %d; %s,%s" % (sat_indx,start_str,end_str)
        
        return out_string

    def __repr__(self):
        if not self.scheduled_dv:
            return  '('+self.get_route_string()+')'
        else:
            return  '('+self.get_route_string()+'; sched dv: %.0f/%.0f Mb'%( self.scheduled_dv, self.data_vol)+')'

    def validate_route (self):
        """ validates timing and ordering of route
        
        [description]
        :raises: Exception, Exception, Exception
        """

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
                raise Exception(string)
            if not wind.start >= last_time or not wind.end >= last_time:
                string ='routing_objects.py: Found a bad start time at window indx %d in route. Route string: %s'%( windex, self.get_route_string())
                raise Exception( string)

            #  note that we manually trace the satellite index through cross-link window here. This is a bit redundant with the functionality of window_start_sats,  but adds a little bit more of a warm, happy, comfortable feeling in the area checking
            if type (wind) is XlnkWindow:
                curr_sat_indx = next_sat_indx
                next_sat_indx = wind.xsat_indx if not wind.xsat_indx == curr_sat_indx else wind.sat_indx

            last_time = wind.end



class LinkInfo(object):
    """docstring fos  LinkInfo"""
    def __init__(self,data_routes=[],total_data_vol=0,used_data_vol=0):
        self.data_routes = data_routes
        self.total_data_vol = total_data_vol
        self.used_data_vol = used_data_vol

    def __str__( self):
        return  "routes: "+str(self.data_routes) + " ; dv %.0f/%.0f Mb" % ( self.used_data_vol, self.total_data_vol)


