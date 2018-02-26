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

    def  sort_windows(self):
        self.route.sort(key=lambda x: x.start)

        #  do some basic validation that the windows are in order
        last_wind_end = const.UNASSIGNED_DT_NEG_INF
        for wind in self.route:
            if wind.start <last_wind_end  or wind.end  < last_wind_end:
                raise Exception ('routing_objects.py:  windows are out of order')
            last_wind_end = wind.end

    def  print_route( self,  time_base= None):
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
        
        print (out_string)

class LinkInfo(object):
    """docstring fos  LinkInfo"""
    def __init__(self,data_routes=[],total_data_vol=0,used_data_vol=0):
        self.data_routes = data_routes
        self.total_data_vol = total_data_vol
        self.used_data_vol = used_data_vol

    def __str__( self):
        return  "routes: "+str(self.data_routes) + " ; dv %.0f/%.0f Mb" % ( self.used_data_vol, self.total_data_vol)


