from datetime import timedelta

from circinus_tools  import  constants as const

class ActivityWindow(object):
    """ specifies an activity that occurs from a start to an end time
    
    note that the hash function for this class uses the window_ID attribute. This should be globally unique across all  activity window instances and subclass instances created in the simulation
    """

    def __init__(self, start, end, window_ID):
        '''
        Creates an activity window

        :param datetime start: start time of the window
        :param datetime end: end time of the window
        :param int window_ID:  unique window ID used for hashing and comparing windows
        '''

        self.start = start
        self.end = end
        self.window_ID = window_ID
        self.data_vol = const.UNASSIGNED
        self.scheduled_data_vol = const.UNASSIGNED
        self.remaining_data_vol = const.UNASSIGNED

    def __hash__(self):
        return self.window_ID

    def __eq__(self, other):
        return self.window_ID ==  other.window_ID

    def update_duration_from_scheduled_dv( self):
        """ update duration based on schedule data volume
        
        updates the schedule duration for the window based upon the assumption that the data volume scheduled for the window is able to be transferred at an average data rate. Updated window times are based off of the center time of the window.
        """

        center_time = self.start + (self.end - self.start)/2

        average_data_rate = self.data_vol/( self.end- self.start).total_seconds()

        scheduled_time_s = self.scheduled_data_vol/average_data_rate

        self.start = center_time - timedelta ( seconds = scheduled_time_s/2)
        self.end = center_time + timedelta ( seconds = scheduled_time_s/2)


    def print_self(self):
        print('ActivityWindow')
        print('start: ' + str(self.start))
        print('end: ' + str(self.end))
        print('......')

    def combine_with_window(self,other_act):
        '''
        Note: deprecated

        Combine this window with another one. The combined window object is stored in self. Overlap existence is assumed -compareWindows should be called first to see if there's overlap or not

        :param other_act: other window to combine into self
        :return: nothing (self stores combined window)
        '''

        if self.start < other_act.start:
            if self.end < other_act.end:
                self.end = other_act.end
            else:
                print('ActivityWindow.py: warning, combining non-overlapping windows')

        if self.start > other_act.start:
            self.start = other_act.start

            if self.start > other_act.end:
                print('ActivityWindow.py: warning, combining non-overlapping windows')

