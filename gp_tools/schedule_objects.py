from copy import copy, deepcopy
from datetime import timedelta

class Dancecard(object):
    def __init__(self, dancecard_start, dancecard_end, tstep_sec):
        '''
        An entity schedule object maintains a schedule data structure that keeps track of the activities that an entity (e.g. satellite, ground station) has committed to perform, and the start and stop times for those acts

        Note that N is equal to the total duration of the dance card divided by the time step in seconds. This is the number of TIMESTEPS, that is, time deltas. the number of distinct TIMEPOINTS in the dance card is time steps +1

        Props to Emily Clements for the idea for the name of the data struct

        Dancecard structure:
        |  i=0    |  i=1    |  i=2    |  i=3    ...   |  i=N    |
        | tstep 0 | tstep 1 | tstep 2 | tstep 3 ...   | tstep N |
        ^         ^         ^         ^               ^         ^ 
        t_pnt0    t_pnt1    t_pnt2    t_pnt3          t_pnt4    t_pnt(N+1)
        ^                             ^
        base_time                *    t    *

        :param scheduling_start: the start of the scheduling window (generally same as scenario)
        :param scheduling_end: the end of the scheduling window (generally same as scenario)
        :param tstep_sec: time step in seconds in overall scenario
        '''

        self.total_duration = (dancecard_end - dancecard_start).total_seconds()
        # each index in dancecard represents a chunk of time, e.g. [1] is from start_time+tstep_sec*1 to start_time+tstep_sec*2
        num_timesteps = int(self.total_duration / tstep_sec)

        # this "dancecard" stores the objects for a given index
        self.dancecard = [[] for i in range(num_timesteps)]

        self.dancecard_start = dancecard_start
        self.dancecard_end = dancecard_end
        self.tstep_sec = tstep_sec
        self.tstep_td = timedelta(seconds=tstep_sec)
        self.num_timesteps = num_timesteps
        self.num_timepoints = num_timesteps+1

    def get_timestep_indices ( self):
        """ return the list of time step indices used for this dance card
        
        returns a list of indices that are valid for every given time step within this dance card. this is useful for doing other operations that need to refer to the same indices as this dance card
        """
        return  range (self.num_timesteps)

    def get_timepoint_indices ( self):
        """ return the list of time point indices used for this dance card
        
        returns a list of indices that are valid for every given time point within this dance card. this is useful for doing other operations that need to refer to the same indices as this dance card
        """
        # TODO: this method should probably be a generator

        return  range (self.num_timepoints)

    def get_timepoint_values ( self,units='seconds',time_option ='relative_to_start'):
        """ return the list of indices used for this dance card
        
        returns a list of indices that are valid for every given time step within this dance card. this is useful for doing other operations that need to refer to the same indices as this dance card
        """

        # TODO: this method should probably be a generator

        if time_option == 'relative_to_start':
            pass
        # TODO: should implement this option - helpful for case where you want to get datetimes back
        elif time_option == 'absolute':
            raise NotImplementedError
        else:
            raise NotImplementedError

        t_vals_s = [self.tstep_sec*tp_indx for tp_indx in range(self.num_timepoints)]

        if units == 'seconds':
            return t_vals_s
        elif units == 'minutes':
            # convert like this to minimize the chance of error buildup
            return [t_val/60.0 for t_val in t_vals_s]


    def get_objects_pre_timepoint_indx(self,tp_indx):
        """ get objects for timestep preceding timepoint index tp_indx
        
        Gets the objects stored in the dance card during the timestep immediately preceding the timepoint represented by tp_indx
        :param tp_indx:  index for the timepoint of interest
        :type tp_indx: int
        :returns: objects at index
        :rtype: {list}
        """

        #  the index of the preceding time step is always tp_indx - 1 
        if tp_indx == 0:
            raise ValueError('Cannot get proceeding timestep for timepoint index 0')

        return self.dancecard[tp_indx-1]

    def get_objects_at_timestep(self,ts_indx):
        """ get objects stored in the dance card at timestep index
        
        Gets the objects stored in the dance card during the timestep at index ts_indx. 
        :param ts_indx:  index for the timestep of interest
        :type ts_indx: int
        :returns: objects at index
        :rtype: {list}
        """

        return self.dancecard[ts_indx]

    @staticmethod
    def get_post_index_dancecard(toi, base_time, tstep_sec):
        '''
        Get the index of the slice of time (of duration tstep_sec) immediately after input time of interest (toi), assuming a dancecard of discrete timesteps. Will first round toi to the nearest timepoint based on tstep_sec.

        e.g.
        Dancecard structure:
        |  i=0    |  i=1    |  i=2    |  i=3    ...   |  i=N    |
        | tstep 0 | tstep 1 | tstep 2 | tstep 3 ...   | tstep N |
                            ^
                       *   toi   * (toi anywhere btwn *) --> returns 2
        ^
        base_time

        :param datetime toi: time of interest for which we want to find the immediately-following chunk of time of an assumed dancecard
        :param datetime base_time: time right before first element/time slice of dancecard
        :param float tstep_sec: the size of timestep in the dancecard
        :return: the index of the timeslice AFTER toi. Normally you should subtract 1 from this if you're looking for an activity end index
        '''

        toi_sec = (toi - base_time).total_seconds()
        indx = int(round(toi_sec / tstep_sec))

        return indx

    def get_post_index_by_time(self, toi):
        '''
        Get the index of the slice of time (of duration tstep_sec) immediately after input time of interest (toi), using the base time for this dancecard. Will first round toi to the nearest timepoint based on tstep_sec.

        :param toi: time of interest for which we want to find the immediately-following chunk of time in this dancecard
        :return: the index of the timeslice AFTER toi. Normally you should subtract 1 from this if you're looking for an activity end index
        '''
        return Dancecard.get_post_index_dancecard(toi, self.dancecard_start, self.tstep_sec)

    def get_pre_time_by_index(self, indxoi):
        '''
        Get the absolute time immediately before a corresponding timestep index, using the base time for this dancecard

        :param int indxoi: index of interest
        :return: the absolute time
        '''

        return self.dancecard_start + timedelta(seconds=indxoi*self.tstep_sec)

    def add_winds_to_dancecard(self, winds):
        '''
        Add a set of windows to the dancecard according to their start and stop times

        :param winds: list of windows to add. Can be a list with a single element, of course
        :return:
        '''
        dancecard_last_indx = len(self.dancecard) - 1

        for wind in winds:
            act_start_indx = Dancecard.get_post_index_dancecard(wind.start, self.dancecard_start, self.tstep_sec)
            act_end_indx = Dancecard.get_post_index_dancecard(wind.end, self.dancecard_start, self.tstep_sec) - 1 # -1 because last time slice is before wind.end

            act_start_indx = max(0, act_start_indx)
            act_end_indx = min(dancecard_last_indx, act_end_indx)

            for indx in range(act_start_indx, act_end_indx + 1): # go all the way up to end indx
                self.dancecard[indx].append(wind)  # add object to schedule. Note this is not a deepcopy!

    def remove_winds_from_dancecard(self, winds, unmodified_yes = False):
        '''
        Remove a set of windows from the dancecard according to their start and stop times. Remove only the objects corresponding to the winds from the dancecard

        :param winds: list of windows to remove. Can be a list with a single element, of course
        :param bool unmodified_yes: set to true if you want to remove everything from the initially created start to end time (unmodded by search process)
        :return:
        '''
        dancecard_last_indx = len(self.dancecard) - 1

        for wind in winds:
            wind_start = wind.start
            wind_end = wind.end
            if unmodified_yes:
                wind_start = wind.unmodified_start
                wind_end = wind.unmodified_end

            act_start_indx = Dancecard.getPostIndexDancecard(wind_start, self.dancecard_start, self.tstep_sec)
            act_end_indx = Dancecard.getPostIndexDancecard(wind_end, self.dancecard_start, self.tstep_sec) - 1 # -1 because last time slice is before wind.end

            act_start_indx = max(0, act_start_indx)
            act_end_indx = min(dancecard_last_indx, act_end_indx)

            for indx in xrange(act_start_indx, act_end_indx + 1): # go all the way up to end indx

                # todo: this is super hacky, shouldn't use try...catch as a shortcut for normal operations
                try:
                    self.dancecard[indx].remove(wind)  # remove object
                except ValueError:
                    pass  #already been removed

    def remove_wind_durations_from_dancecard(self, winds, unmodified_yes = False):
        '''
        Remove a set of windows as well as all potentially overlapping windows from the dancecard according to their start and stop times. Removes not only the objects from winds, but also everything else in the dancecard for the durations of the winds.

        Note that this function is dangerous to call on a dancecard if you don't know what you're doing. For example, crosslink objects appear in the dancecards for both crosslinking sats. Could leave one sat thinking there's still xlnk availability when there's not

        :param winds: list of windows to remove. Can be a list with a single element, of course
        :param bool unmodified_yes: set to true if you want to remove everything from the initially created start to end time (unmodded by search process)
        :return:
        '''

        dancecard_last_indx = len(self.dancecard) - 1

        for wind in winds:
            wind_start = wind.start
            wind_end = wind.end
            if unmodified_yes:
                wind_start = wind.unmodified_start
                wind_end = wind.unmodified_end

            act_start_indx = Dancecard.get_post_index_dancecard(wind_start, self.dancecard_start, self.tstep_sec)
            act_end_indx = Dancecard.get_post_index_dancecard(wind_end, self.dancecard_start, self.tstep_sec) - 1 # -1 because last time slice is before wind.end

            act_start_indx = max(0, act_start_indx)
            act_end_indx = min(dancecard_last_indx, act_end_indx)

            for indx in range(act_start_indx, act_end_indx + 1): # go all the way up to end indx
                self.dancecard[indx] = [] # clear all objects from the dancecard here.