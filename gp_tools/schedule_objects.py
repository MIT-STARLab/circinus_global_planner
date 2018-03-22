from copy import copy, deepcopy
from datetime import timedelta

class Dancecard(object):
    def __init__(self, dancecard_start_dt, dancecard_end_dt, tstep_sec, item_init=list):
        """ Maintains a time series of objects for use in scheduling problems

        Note that N is equal to the total duration of the dance card divided by 
        the time step in seconds. This is the number of TIMESTEPS, that is, 
        time deltas. the number of distinct TIMEPOINTS in the dance card is time steps +1

        Props to Emily Clements for the idea for the name of the data struct

        Dancecard structure:
        |  i=0    |  i=1    |  i=2    |  i=3    ...   |  i=N    |
        | tstep 0 | tstep 1 | tstep 2 | tstep 3 ...   | tstep N |
        ^         ^         ^         ^               ^         ^ 
        tpnt0     tpnt1     tpnt2     tpnt3           tpnt4     tpnt(N+1)
        ^                             ^
        base_time                *    t    *

        :param dancecard_start_dt:  start time for the dance card
        :type dancecard_start_dt: datetime
        :param dancecard_end_dt: end time for the dance card
        :type dancecard_end_dt: datetime
        :param tstep_sec:  time step in seconds
        :type tstep_sec:  float
        :param item_init: how to initialize indices in card, defaults to list
        :type item_init: {list,None}, optional
        """

        self.total_duration = (dancecard_end_dt - dancecard_start_dt).total_seconds()
        # each index in dancecard represents a chunk of time, e.g. [1] is from start_time+tstep_sec*1 to start_time+tstep_sec*2
        num_timesteps = int(self.total_duration / tstep_sec)

        if item_init == list:
            # this "dancecard" stores a list of objects for a given index
            self.dancecard = [[] for i in range(num_timesteps)]
        elif item_init == None:
            # this "dancecard" stores for a given index
            self.dancecard = [None for i in range(num_timesteps)]
        else:
            raise NotImplementedError

        self.dancecard_start_dt = dancecard_start_dt
        self.dancecard_end_dt = dance_dtcard_end
        self.tstep_sec = tstep_sec
        self.tstep_td = timedelta(seconds=tstep_sec)
        self.num_timesteps = num_timesteps
        self.num_timepoints = num_timesteps+1

    def get_ts_indcs ( self):
        """ return time step indices in this dance card
        
        returns a list of time step indices for this dance card. These can be used directly for accessing the objects stored in a time step. each time step falls between two time points on either side
        :returns:  list of all time indcs in dance card
        :rtype: {[list of int]}
        :raises: NotImplementedError
        """

        return  range (self.num_timesteps)

    def get_tp_indcs ( self):
        """ return time point indices in this dance card
        
        returns a list of time point indices for this dance card. these are indices that correspond exactly to the time point values in the dance card. they cannot be used directly for accessing the objects stored in a time step, but they tell you the bounds of each time step, and can be used for accessing the objects either before or after a given time point. This is useful for doing other operations that need to refer to the same times as this dance card
        :returns:  list of all time indcs in dance card
        :rtype: {[list of int]}
        :raises: NotImplementedError
        """

        # TODO: this method should probably be a generator

        return  range (self.num_timepoints)

    def get_tp_values ( self,out_units='seconds'):
        """ return time points in this dance card
        
        returns a list of time points for this dance card. this is useful for doing other operations that need to refer to the same times as this dance card
        :returns:  list of all time points in dance card
        :rtype: {[list, element type specified by out_units]}
        :raises: NotImplementedError
        """

        # TODO: this method should probably be a generator


        t_vals_s = [self.tstep_sec*tp_indx for tp_indx in range(self.num_timepoints)]

        if units == 'seconds':
            return t_vals_s
        elif units == 'minutes':
            # convert like this to minimize the chance of error buildup
            return [t_val/60.0 for t_val in t_vals_s]
        else:
            raise NotImplementedError

    def get_tp_from_tp_indx(self,tp_indx,out_units='seconds'):
        """get time point value from time point index
        
        [description]
        :param tp_indx:  time point index
        :type tp_indx: int
        :param out_units:  type for return value, defaults to 'seconds'
        :type out_units: str, optional
        :returns:  time point value
        :rtype: {[ type specified by out_units]}
        :raises: NotImplementedError
        """

        if units == 'seconds':
            return self.tstep_sec*tp_indx 
        elif units == 'minutes':
            return self.tstep_sec*tp_indx/60.0
        elif units == 'datetime':
            return self.dancecard_start_dt + timedelta(seconds= self.tstep_sec*tp_indx)
        else:
            raise NotImplementedError

    def get_tp_indx_pre_t(self,t,in_units='datetime'):
        """ get closest time point index preceding time t
        
        [description]
        :param t:  a time within the start and end of the dance card
        :type t: [specified by in_units]
        :param in_units: type for t, defaults to 'datetime'
        :type in_units: str, optional
        :returns:  time point index
        :rtype: {int}
        :raises: NotImplementedError
        """

        if units == 'datetime':
            return floor((t - self.dancecard_start_dt).total_seconds() / self.tstep_sec)
        else:
            raise NotImplementedError

    def get_objects_pre_tp_indx(self,tp_indx):
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

    def get_objects_at_ts_indx(self,ts_indx):
        """ get objects stored in the dance card at timestep index
        
        Gets the objects stored in the dance card during the timestep at index ts_indx. 
        :param ts_indx:  index for the timestep of interest
        :type ts_indx: int
        :returns: objects at index
        :rtype: {list}
        """

        return self.dancecard[ts_indx]

    @staticmethod
    def get_post_ts_indx(toi, base_time, tstep_sec):
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

    def get_post_ts_indx_from_t(self, t, in_units='datetime'):
        """get index of time step following the time point closest to t
        
        Get the index of the slice of time (of duration tstep_sec) immediately after input time of interest (t), using the base time for this dancecard. Will first round t to the nearest timepoint based on tstep_sec.

        :param t:  the time of interest
        :type t: [specified by in_units]
        :param in_units: type for t, defaults to 'datetime'
        :type in_units: str, optional
        :returns:  time step index
        :rtype: {int}
        :raises: NotImplementedError
        """

        if in_units == 'datetime':
            return Dancecard.get_post_ts_indx(t, self.dancecard_start_dt, self.tstep_sec)
        else:
            raise NotImplementedError

    def get_pre_tp_from_ts_indx(self, ts_indx,out_units='datetime'):
        """ get time point preceding time step index
        
        [description]
        :param ts_indx:  time step index
        :type ts_indx: int
        :param out_units:  type for return value, defaults to 'datetime'
        :type out_units: str, optional
        :returns:  time point
        :rtype: {[specified by out_units]}
        :raises: NotImplementedError
        """

        if out_units == 'datetime':
            return self.dancecard_start_dt + timedelta(seconds=ts_indx*self.tstep_sec)
        else:
            raise NotImplementedError

    def add_winds_to_dancecard(self, winds):
        '''
        Add a set of windows to the dancecard according to their start and stop times

        :param winds: list of windows to add. Can be a list with a single element, of course
        :return:
        '''
        dancecard_last_indx = len(self.dancecard) - 1

        for wind in winds:
            act_start_indx = Dancecard.get_post_ts_indx(wind.start, self.dancecard_start_dt, self.tstep_sec)
            act_end_indx = Dancecard.get_post_ts_indx(wind.end, self.dancecard_start_dt, self.tstep_sec) - 1 # -1 because last time slice is before wind.end

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

            act_start_indx = Dancecard.get_post_ts_indx(wind_start, self.dancecard_start_dt, self.tstep_sec)
            act_end_indx = Dancecard.get_post_ts_indx(wind_end, self.dancecard_start_dt, self.tstep_sec) - 1 # -1 because last time slice is before wind.end

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

            act_start_indx = Dancecard.get_post_ts_indx(wind_start, self.dancecard_start_dt, self.tstep_sec)
            act_end_indx = Dancecard.get_post_ts_indx(wind_end, self.dancecard_start_dt, self.tstep_sec) - 1 # -1 because last time slice is before wind.end

            act_start_indx = max(0, act_start_indx)
            act_end_indx = min(dancecard_last_indx, act_end_indx)

            for indx in range(act_start_indx, act_end_indx + 1): # go all the way up to end indx
                self.dancecard[indx] = [] # clear all objects from the dancecard here.