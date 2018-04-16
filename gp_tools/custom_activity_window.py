# Code for dealing with the types of activity windows used in scheduling

# @author: Kit Kennedy

# a note on indices and IDs: indx is the internal representation of an object's ID -  it encodes location of the objects in a number of data structures that hold data for the full set of objects. ID is meant to encode a human readable ID for an object, which may or may not be the same as the internally-used indx. Processing code should be written to operate on indx, not on ID,  unless it's dealing with input or output.

from datetime import datetime, timedelta

from numpy import mean as np_mean

from circinus_tools  import time_tools as tt
from circinus_tools  import  constants as const
from circinus_tools.activity_window import ActivityWindow

DATE_STRING_FORMAT = 'short'
# DATE_STRING_FORMAT = 'iso'

def short_date_string(dt):
    return dt.strftime("%H:%M:%S")

def date_string(dt):
    if DATE_STRING_FORMAT == 'iso':
        return dt.isoformat()
    if DATE_STRING_FORMAT == 'short':
        return  short_date_string(dt)

class ObsWindow(ActivityWindow):
    def __init__(self, window_ID, sat_indx, target_IDs, sat_target_indx, is_urgent, start, end):
        '''
        An observation window. Can represent a window during which an activity can happen, or the actual activity itself

        :param int window_ID: a unique ID for this obs window, for hashing and comparing windows
        :param int sat_indx: index of the satellite
        :param list target_IDs: the set of target IDs that this observation is looking at
        :param int sat_target_indx: the index of this observation in the list of sat-target observations from gp_input_obs.mat
        :param bool is_urgent: set true if the obs is urgent, and it should be routed before all others
        :param datetime start: start time of the window
        :param datetime end: end time of the window
        '''

        self.sat_indx = sat_indx
        self.target_IDs = target_IDs
        self.sat_target_indx = sat_target_indx
        self.data_pkts = []
        self.is_urgent = is_urgent
        self.unmodified_data_vol = const.UNASSIGNED
        self.collected_data_vol = 0
        super(ObsWindow, self).__init__(start, end, window_ID)

    def print_self(self):
        print('ObsWindow')
        print('sat_indx: ' + str(self.sat_indx))
        print('target_IDs: ' + str(self.target_IDs))
        print('window_ID: ' + str(self.window_ID))
        print('start: ' + str(self.start))
        print('end: ' + str(self.end))
        print('......')

    def combine_with_window(self, other_obs):
        for target_ID in other_obs.target_IDs:
            self.target_IDs.append(target_ID)

        super(ObsWindow, self).combine_with_window(other_obs)

    def calc_data_vol(self,pl_data_rate):
        self.data_vol = (self.end- self.start).total_seconds()*pl_data_rate
        self.remaining_data_vol = self.data_vol
        self.unmodified_data_vol = self.data_vol

    def __str__(self):
        return  "(ObsWindow id %d sat %d dv %f targs %s %s,%s)" % ( self.window_ID, self.sat_indx,  self.data_vol,str(self.target_IDs),date_string(self.start),date_string(self.end))

    def __repr__(self):
        return  "(ObsWindow id %d sat %d dv %f targs %s %s,%s)" % (self.window_ID,self.sat_indx,  self.data_vol,str(self.target_IDs),date_string(self.start),date_string(self.end))


class CommWindow(ActivityWindow):
    def __init__(self, start, end,window_ID):
        super(CommWindow, self).__init__(start, end,window_ID)

    def set_data_vol(self,rates_mat,rates_mat_dv_indx=1):
        """
        Calculates the total data volume that can be sent over this link. Uses average data rate to determine data volume. Depending on how much the input data rates matrix is decimated, this could lead to over or underestimates of data volume.

        :return:
        """

        # Note: float[num_timepoints][2] rates_mat: matrix of datarates at each time during the pass. First column is time in MJD, and second column is data rate from sat to xsat in Mbps, third is rate from xsat to sat.

        start_mjd = tt.datetime2mjd(self.start)-5/86400.0  # add 5 secs padding to evade any precision problems
        end_mjd = tt.datetime2mjd(self.end)+5/86400.0  # add 5 secs padding to evade any precision problems

        #  this is fixed in the structure of the data rates output file
        rates_mat_tp_indx = 0;

        data_rates = []
        for i in range(len(rates_mat)):
            # if point i is within window -  this should take care of any indexing issues
            if rates_mat[i][rates_mat_tp_indx] >= start_mjd and rates_mat[i][rates_mat_tp_indx] <= end_mjd:
                data_rates.append(rates_mat[i][rates_mat_dv_indx])

        try:
            #  take the average of all the data rates we saw and multiply by the duration of the window to get data volume
            self.data_vol = np_mean(data_rates) * (self.end - self.start).total_seconds()
        except RuntimeWarning as e:
            raise RuntimeWarning('Trouble determining average data rate. Probable no time points were found within start and end of window. Ensure that you are not overly decimating data rate calculations in data rates input file (window: %s, exception seen: %s)'%(self,str(e)))


class DlnkWindow(CommWindow):
    def __init__(self, window_ID, sat_indx, gs_indx, sat_gs_indx, start, end):
        '''
        A downlink window. Can represent a window during which an activity can happen, or the actual activity itself

        :param int window_ID: a unique ID for this window, for hashing and comparing windows
        :param int sat_indx: index of the satellite
        :param int gs_indx: ground station indx for this downlink
        :param int sat_gs_indx: the index of this dlnk in the list of sat-gs downlinks from gp_input_gslink.mat
        :param datetime start: start time of the window
        :param datetime end: end time of the window
        '''

        self.sat_indx = sat_indx
        self.gs_indx = gs_indx
        self.sat_gs_indx = sat_gs_indx
        self.routed_data_vol = 0

        super(DlnkWindow, self).__init__(start, end, window_ID)

    # THIS IS A DIRTY HACK for dealing with a legacy pick pickle file that uses gs_ID instead of gs_indx. 
    # TODO: remove this after  initial dev
    def __getattr__(self, attr):
        if attr == 'gs_indx' and 'gs_indx' not in dir(self):
            return self.gs_ID
        else:
            super().__getattr__(attr)            

    def print_self(self,  print_data_vol = True):
        print('DlnkWindow')
        print('sat_indx: ' + str(self.sat_indx))
        print('gs_indx: ' + str(self.gs_indx))
        print('sat_gs_indx: ' + str(self.sat_gs_indx))
        print('window_ID: ' + str(self.window_ID))
        if print_data_vol:  print('data_vol: ' + str(self.data_vol))
        print('start: ' + str(self.start))
        print('end: ' + str(self.end))
        print('......')

    def __str__(self):
        return  "(DlnkWindow id %d sat %d dv %f gs %d %s,%s)" % (self.window_ID,self.sat_indx,  self.data_vol, self.gs_indx,date_string(self.start),date_string(self.end))

    def __repr__(self):
        return  "(DlnkWindow id %d sat %d dv %f gs %d %s,%s)" % (self.window_ID,self.sat_indx,  self.data_vol, self.gs_indx,date_string(self.start),date_string(self.end))

class XlnkWindow(CommWindow):
    def __init__(self, window_ID, sat_indx, xsat_indx, sat_xsat_indx, start, end, symmetric=True, tx_sat=None):
        '''
        A downlink window. Can represent a window during which an activity can happen, or the actual activity itself

        A cross-link window can be marked as symmetric as a space-saving mechanism. "symmetric" means that both satellites are performing the cross-link window with the same start and end and at the same data rate ( given that the satellites have the same communication technology on both sides, both are equally able to transmit and receive, and there are no unidirectional effects in the link, this will be the case). This saves space by only using one object to represent the cross-link.  it's almost like I thought about this really hard in advance and realized that this was a good way to go - nah, I added in directionality later.

        If the window is not symmetric, then one of the satellite should be marked as the transmitting satellite. note that the sat_indx and xsat_indx fields should not be used to infer the direction of the cross-link.

        :param int window_ID: a unique ID for this window, for hashing and comparing windows
        :param int sat_indx: index of the satellite on "side A" of the link
        :param int xsat_indx: index of the satellite on "side B" of the link (should be > sat_indx)
        :param int sat_xsat_indx: the serial index of this xlnk in the list of sat-xsat links from input struct
        :param datetime start: start time of the window
        :param datetime end: end time of the window
        :param bool symmetric:  true if both satellites are transmitting at the same data rate during the window.
        :param bool tx_sat: the single satellite that is transmitting, if the window is not symmetric
        '''

        self.sat_indx = sat_indx
        self.xsat_indx = xsat_indx
        self.sat_xsat_indx = sat_xsat_indx
        self.symmetric = symmetric
        #  note that tx_sat is NOT meant to be used to keep track of which satellite was decided on as the tx-er for a symmetric window (e.g. after route selection). Should keep track of that in an external data structure (e.g. in a data route object), otherwise there's a lot of potential for problems to crop up from the fact that this cross-link window object may appear in multiple routes
        self.tx_sat = tx_sat

        if symmetric and tx_sat is not None:
            raise RuntimeError('Cross-link window should not be both symmetric and have a transmitting satellite specified')
        if not symmetric and tx_sat is None:
            raise RuntimeError('Cross-link window should either be marked as symmetric or have a transmitting satellite specified')

        super(XlnkWindow, self).__init__(start, end, window_ID)

    def print_self(self,  print_data_vol = True):
        print('XlnkWindow')
        print('sat_indx: ' + str(self.sat_indx))
        print('xsat_indx: ' + str(self.xsat_indx))
        print('sat_xsat_indx: ' + str(self.sat_xsat_indx))
        print('window_ID: ' + str(self.window_ID))
        if print_data_vol:  print('data_vol: ' + str(self.data_vol))
        print('start: ' + str(self.start))
        print('end: ' + str(self.end))
        print('......')

    def __str__(self):
        return  self.__repr__()

    def __repr__(self):
        if self.symmetric:
            return  "(XlnkWindow id %d sym. sats %d<->%d dv %f %s,%s)" % (self.window_ID,self.sat_indx, self.xsat_indx,self.data_vol, date_string(self.start),date_string(self.end))
        else:
            return  "(XlnkWindow id %d uni. sats %d->%d dv %f %s,%s)" % (self.window_ID,self.tx_sat, self.get_xlnk_partner(self.tx_sat),self.data_vol, date_string(self.start),date_string(self.end))

    def get_xlnk_partner(self,sat_indx):
        """return indx of the cross-linked partner
        
        This method tells us, from sat_indx's perspective, who its cross-link partner is. Does not take directionality of window into account
        :param sat_indx: a satellite index contained in either self.sat_indx or self.xsat_indx
        :type sat_indx: int
        :returns: the index of the crosslink partner satellite
        :rtype: {int}
        """

        xsat_indx= self.xsat_indx  if sat_indx == self.sat_indx else self.sat_indx
        return xsat_indx

    def is_rx(self,sat_indx):
        """ check if satellite is receiving during this cross-link
        
        satellite is always a receiver during a symmetric cross-link, and is a receiver if it's not the tx satellite in a nonsymmetric cross-link
        :param sat_indx:  satellite index
        :type sat_indx: int
        :returns:  true if it's receiving
        :rtype: {bool}
        """
        if self.symmetric:
            return True
        elif sat_indx == self.tx_sat:
            return False
        else:
            return True

class UrgentWindow(ActivityWindow):
    def __init__(self, target_ID, start, end):
        '''
        Window of time during which an observation target is urgent. Convenience class.

        :param int target_ID: target ID of target to make urgent during this window
        :param datetime start: start time of the window
        :param datetime end: end time of the window
        '''

        self.target_ID = target_ID
        super(UrgentWindow, self).__init__(start, end)


class EclipseWindow(ActivityWindow):
    def __init__(self, window_ID, start, end):
        '''
        An eclipse window. Meant to represent when a satellite is in eclipse and can't see the sun

        :param int window_ID: a unique ID for this obs window, for hashing and comparing windows
        :param datetime start: start time of the window
        :param datetime end: end time of the window
        '''

        super(EclipseWindow, self).__init__(start, end, window_ID)

    def __str__(self):
        return  "(EclipseWindow id %d %s,%s)" % ( self.window_ID,date_string(self.start),date_string(self.end))

    def __repr__(self):
        return  "(EclipseWindow id %d %s,%s)" % ( self.window_ID,date_string(self.start),date_string(self.end))
