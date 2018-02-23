# Code for dealing with the types of activity windows used in scheduling

# @author: Kit Kennedy

from datetime import datetime, timedelta

from circinus_tools  import time_tools as tt
from circinus_tools  import  constants as const
from circinus_tools.activity_window import ActivityWindow


class ObsWindow(ActivityWindow):
    def __init__(self, window_ID, sat_indx, target_IDs, sat_target_indx, is_urgent, start, end):
        '''
        An observation window. Can represent a window during which an activity can happen, or the actual activity itself

        :param int window_ID: a unique ID for this obs window, for ease of tracking
        :param int sat_indx: index of the satellite
        :param list target_IDs: the set of target IDs that this observation is looking at
        :param int sat_target_indx: the index of this observation in the list of sat-target observations from gp_input_obs.mat
        :param bool is_urgent: set true if the obs is urgent, and it should be routed before all others
        :param datetime start: start time of the window
        :param datetime end: end time of the window
        '''

        self.sat_indx = sat_indx
        self.window_ID = window_ID
        self.target_IDs = target_IDs
        self.sat_target_indx = sat_target_indx
        self.data_pkts = []
        self.is_urgent = is_urgent
        self.data_vol = const.UNASSIGNED
        self.remaining_data_vol = const.UNASSIGNED
        self.unmodified_data_vol = const.UNASSIGNED
        self.collected_data_vol = 0
        super(ObsWindow, self).__init__(start, end)

    def print_self(self):
        print('ObsWindow')
        print('sat_indx: ' + str(self.sat_indx))
        print('target_IDs: ' + str(self.target_IDs))
        print('window_ID: ' + str(self.window_ID))
        print('start: ' + str(self.start))
        print('end: ' + str(self.end))
        print('duration: ' + str(self.duration))
        print('......')

    def combine_with_window(self, other_obs):
        for target_ID in other_obs.target_IDs:
            self.target_IDs.append(target_ID)

        super(ObsWindow, self).combine_with_window(other_obs)

    def calc_data_vol(self,pl_data_rate):
        self.data_vol = self.duration.total_seconds()*pl_data_rate
        self.remaining_data_vol = self.data_vol
        self.unmodified_data_vol = self.data_vol

class CommWindow(ActivityWindow):
    def __init__(self, start, end):
        self.data_vol = const.UNASSIGNED
        self.remaining_data_vol = const.UNASSIGNED
        self.data_pkts = []
        # store initial start and end so we have them after start and end themselves get modified
        self.unmodified_start = start
        self.unmodified_end = end
        self.unmodified_data_vol = const.UNASSIGNED
        super(CommWindow, self).__init__(start, end)

    def set_data_vol_and_refresh_times(self):
        """
        Calculates the total data volume that can be sent over this link. Resets the start and stop times based on the non-zero rates in rates_mat. If the rates mat has all zeros as the data rates, the data volume is set to 0. Note that the original duration of the window (self) will be collapsed to hug the non-zero data rate points in rates_mat. If the rate values in rates mat are undersampled such that they don't exactly line up with start and stop times of original window, the window will be shortened a little bit.

        :return:
        """

        # Note: float[num_timepoints][2] rates_mat: matrix of datarates at each time during the pass. First column is time in MJD, and second column is data rate in Mbps. It is assumed that the data rate transitions from 0 to non-zero only once, and back from non-zero to 0 only once.

        start_mjd = tt.datetime2mjd(self.start)-5/86400.0  # add 5 secs padding to evade any precision problems
        end_mjd = tt.datetime2mjd(self.end)+5/86400.0  # add 5 secs padding to evade any precision problems

        first_mjd_with_nonzero = const.UNASSIGNED
        last_mjd_with_nonzero = const.UNASSIGNED
        found_first = False
        # found_last = False

        rates_mat = self.rates_mat

        # note that the below code assumes the data rate over the whole interval from time i to time i+1 is the data rate at i.

        total_data_vol = 0
        for i in range(len(rates_mat)):

            # if point i is within window, and has a data rate greater than 0.
            if (rates_mat[i][0] >= start_mjd and rates_mat[i][0] <= end_mjd) and (rates_mat[i][1] > 0):

                if found_first == False:
                    first_mjd_with_nonzero = rates_mat[i][0]
                    found_first = True

                last_mjd_with_nonzero = rates_mat[i][0]

                # if we still have another interval before the end of the window
                if (i<len(rates_mat)-1)  and (rates_mat[i + 1][0] <= end_mjd):
                    time_slice_sec = (rates_mat[i+1][0] - rates_mat[i][0]) * 86400

                    #                   data rate at i     *  time
                    total_data_vol += rates_mat[i][1] * time_slice_sec


        self.data_vol = total_data_vol
        self.remaining_data_vol = self.data_vol
        self.unmodified_data_vol = total_data_vol

        # If the data rates are not all zero and yet the times are const.UNASSIGNED, there's a problem
        if first_mjd_with_nonzero == const.UNASSIGNED and last_mjd_with_nonzero == const.UNASSIGNED and self.data_vol > 0:
            print('problem with recalculating window start and stop time')
            1/0

        # if no problems with const.UNASSIGNED start/end mjd, then update start and end
        else:
            old_start = self.start - timedelta(seconds=5)
            old_end = self.end + timedelta(seconds=5)
            self.start = tt.mjd2datetime(first_mjd_with_nonzero)
            self.end = tt.mjd2datetime(last_mjd_with_nonzero)
            self.refresh_duration()

            # if the window is actually way too short, indicate by making data vol 0
            if self.end - self.start < timedelta(seconds=5):
                self.data_vol = 0


            if self.end == const.UNASSIGNED_DT and self.start != const.UNASSIGNED_DT:
                1 / 0

            # just some sanity checks here...
            if self.start != const.UNASSIGNED_DT and (self.start < old_start or self.start > old_end):
                1/0
            if self.end != const.UNASSIGNED_DT and (self.end < old_start or self.end > old_end):
                1/0

    def recalc_data_vol(self):
        """
        Calculates the total data volume that can be sent over this link, but does not change any of the data fields of object self

        :return: the data volume that can be sent over the link
        """

        # Note: float[num_timepoints][2] rates_mat: matrix of datarates at each time during the pass. First column is time in MJD, and second column is data rate in Mbps. It is assumed that the data rate transitions from 0 to non-zero only once, and back from non-zero to 0 only once.

        start_mjd = tt.datetime2mjd(self.start)-5/86400.0  # add 5 secs padding to evade any precision problems
        end_mjd = tt.datetime2mjd(self.end)+5/86400.0  # add 5 secs padding to evade any precision problems

        rates_mat = self.rates_mat

        # note that the below code assumes the data rate over the whole interval from time i to time i+1 is the data rate at i.

        total_data_vol = 0
        for i in range(len(rates_mat)):

            # if point i is within window, and has a data rate greater than 0.
            if (rates_mat[i][0] >= start_mjd and rates_mat[i][0] <= end_mjd) and (rates_mat[i][1] > 0):

                # if we still have another interval before the end of the window
                if (i<len(rates_mat)-1)  and (rates_mat[i + 1][0] <= end_mjd):
                    time_slice_sec = (rates_mat[i+1][0] - rates_mat[i][0]) * 86400

                    #                   data rate at i     *  time
                    total_data_vol += rates_mat[i][1] * time_slice_sec


        return total_data_vol

class DlnkWindow(CommWindow):
    def __init__(self, window_ID, sat_indx, gs_ID, sat_gs_indx, start, end):
        '''
        A downlink window. Can represent a window during which an activity can happen, or the actual activity itself

        :param int window_ID: a unique ID for this window, for ease of tracking
        :param int sat_indx: index of the satellite
        :param int gs_ID: ground station ID for this downlink
        :param int sat_gs_indx: the index of this dlnk in the list of sat-gs downlinks from gp_input_gslink.mat
        :param datetime start: start time of the window
        :param datetime end: end time of the window
        '''

        self.sat_indx = sat_indx
        self.window_ID = window_ID
        self.gs_ID = gs_ID
        self.sat_gs_indx = sat_gs_indx
        self.routed_data_vol = 0

        super(DlnkWindow, self).__init__(start, end)

    def print_self(self,  print_data_vol = True):
        print('DlnkWindow')
        print('sat_indx: ' + str(self.sat_indx))
        print('gs_ID: ' + str(self.gs_ID))
        print('sat_gs_indx: ' + str(self.sat_gs_indx))
        print('window_ID: ' + str(self.window_ID))
        if print_data_vol:  print('data_vol: ' + str(self.data_vol))
        print('start: ' + str(self.start))
        print('end: ' + str(self.end))
        print('duration: ' + str(self.duration))
        print('......')

class XlnkWindow(CommWindow):
    def __init__(self, window_ID, sat_indx, xsat_indx, sat_xsat_indx, start, end):
        '''
        A downlink window. Can represent a window during which an activity can happen, or the actual activity itself

        :param int window_ID: a unique ID for this window, for ease of tracking
        :param int sat_indx: index of the satellite on "side A" of the link
        :param int xsat_indx: index of the satellite on "side B" of the link (should be > sat_indx)
        :param int sat_xsat_indx: the serial index of this xlnk in the list of sat-xsat links from input struct
        :param datetime start: start time of the window
        :param datetime end: end time of the window
        '''

        self.sat_indx = sat_indx
        self.window_ID = window_ID
        self.xsat_indx = xsat_indx
        self.sat_xsat_indx = sat_xsat_indx
        # this stores the amount of data routed to either sat index during the xlnk. Dictionary keys are integers representing the sat index routed to, and the values are the amount of data routed
        # Using a dictionary here to support eventual bi-directionality of crosslinks
        self.routed_data_vol_to_sat_indx = {}
        self.routed_pkts_to_sat_indx = {}
        self.routed_pkt_ids_to_sat_indx = {}

        super(XlnkWindow, self).__init__(start, end)

    def print_self(self,  print_data_vol = True):
        print('XlnkWindow')
        print('sat_indx: ' + str(self.sat_indx))
        print('xsat_indx: ' + str(self.xsat_indx))
        print('sat_xsat_indx: ' + str(self.sat_xsat_indx))
        print('window_ID: ' + str(self.window_ID))
        if print_data_vol:  print('data_vol: ' + str(self.data_vol))
        print('start: ' + str(self.start))
        print('end: ' + str(self.end))
        print('duration: ' + str(self.duration))
        print('......')

    def get_total_routed_data(self):
        '''
        Magnitude of total data exchanged in this link, regardless of direction

        :return:
        '''

        return sum(self.routed_data_vol_to_sat_indx[key] for key in self.routed_data_vol_to_sat_indx.keys())

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