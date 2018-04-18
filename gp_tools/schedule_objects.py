# Data structure for representing time series data
# 
# Note: maybe this could be done better with Pandas... works for now though
# 
# @author Kit Kennedy

from copy import copy, deepcopy
from datetime import timedelta
from math import floor, ceil

class Dancecard(object):
    def __init__(self, dancecard_start_dt, dancecard_end_dt, tstep_sec, item_init=list,item_type=list,mode='timestep'):
        """ Maintains a time series of objects for use in scheduling problems

        Note that N is equal to the total duration of the dance card divided by 
        the time step in seconds. This is the number of TIMESTEPS, that is, 
        time deltas. the number of distinct TIMEPOINTS in the dance card is time steps +1

        Props to Emily Clements for the idea for the name of the data struct

        There are two different modes. default is "timestep" mode.
        - in "timestep" mode, each time step or delta has a distinct index in the internal dancecard list and can store data. This mode is meant to represent activities that have durations
        - in "timepoint" mode, each absolute time has a distinct index in the internal dancecard list and can store data. This mode is meant to represent states that are valid at a certain point in time.

        DANCECARD STRUCTURE

        Timestep mode:
        |  i=0    |  i=1    |  i=2    |  i=3    ...   |  i=N    |
        | tstep 0 | tstep 1 | tstep 2 | tstep 3 ...   | tstep N |
        ^         ^         ^         ^               ^         ^ 
        tpnt0     tpnt1     tpnt2     tpnt3           tpnt4     tpnt(N+1)
        ^ 
        base_time

        Timepoint mode:
        |  i=0    |  i=1    |  i=2    |  i=3    ...   |  i=N+1   |
        | tpnt 0  | tpnt 1  | tpnt 2  | tpnt 3  ...   | tpnt N+1 |
                  ^         ^         ^               ^          
                  tstep1    tstep2    tstep3          tstep N
            ^
            base_time

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
        num_timepoints = num_timesteps+1

        if mode == 'timestep':
            if item_init == list:
                # this "dancecard" stores a list of objects for a given index
                self.dancecard = [[] for i in range(num_timesteps)]
            elif item_init == None:
                # this "dancecard" stores for a given index
                self.dancecard = [None for i in range(num_timesteps)]
            else:
                raise NotImplementedError
        elif mode == 'timepoint':
            if item_init == list:
                # this "dancecard" stores a list of objects for a given index
                self.dancecard = [[] for i in range(num_timepoints)]
            elif item_init == None:
                # this "dancecard" stores an arbitrary object for a given index
                self.dancecard = [None for i in range(num_timepoints)]
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError

        self.dancecard_start_dt = dancecard_start_dt
        self.dancecard_end_dt = dancecard_end_dt
        self.tstep_sec = tstep_sec
        self.tstep_td = timedelta(seconds=tstep_sec)
        self.num_timesteps = num_timesteps
        self.num_timepoints = num_timepoints
        self.mode = mode
        self.item_type = item_type

    def __setitem__(self, key, value):
        """ setter for internal dancecard by index"""
        self.dancecard[key] = value

    def __getitem__(self, key):
        """ getter for internal dancecard by index"""
        return self.dancecard[key]

    def get_ts_indcs ( self):
        """ return time step indices in this dance card
        
        returns a generator for time step indices for this dance card. These can be used directly for accessing the objects stored in a time step. each time step falls between two time points on either side
        :returns:  generator for all time indcs in dance card
        :rtype: {[int generator]}
        :raises: NotImplementedError
        """

        return (ts_indx for ts_indx in range (self.num_timesteps))

    def get_tp_indcs ( self):
        """ return time point indices in this dance card
        
        returns a generator of time point indices for this dance card. these are indices that correspond exactly to the time point values in the dance card. they cannot be used directly for accessing the objects stored in a time step, but they tell you the bounds of each time step, and can be used for accessing the objects either before or after a given time point. This is useful for doing other operations that need to refer to the same times as this dance card
        :returns:  generator for all time indcs in dance card
        :rtype: {[int generator]}
        :raises: NotImplementedError
        """

        return (tp_indx for tp_indx in range(self.num_timepoints))

    def get_tp_values ( self,out_units='seconds'):
        """ return time points in this dance card
        
        returns a generator of time points for this dance card. this is useful for doing other operations that need to refer to the same times as this dance card
        :returns:  generator for all time points in dance card
        :rtype: {[generator, element type specified by out_units]}
        :raises: NotImplementedError
        """

        t_vals_s = (self.tstep_sec*tp_indx for tp_indx in range(self.num_timepoints))

        if out_units == 'seconds':
            return t_vals_s
        elif out_units == 'minutes':
            # convert like this to minimize the chance of error buildup
            return (t_val/60.0 for t_val in t_vals_s)
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

        if out_units == 'seconds':
            return self.tstep_sec*tp_indx 
        elif out_units == 'minutes':
            return self.tstep_sec*tp_indx/60.0
        elif out_units == 'datetime':
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

        if in_units == 'datetime':
            if t < self.dancecard_start_dt:
                raise ValueError('t (%s) is before dancecard start (%s)'%(t.isoformat(),self.dancecard_start_dt.isoformat()))

            return floor((t - self.dancecard_start_dt).total_seconds() / self.tstep_sec)
        else:
            raise NotImplementedError

    def get_tp_indx_post_t(self,t,in_units='datetime'):
        """ get closest time point index following time t
        
        [description]
        :param t:  a time within the start and end of the dance card
        :type t: [specified by in_units]
        :param in_units: type for t, defaults to 'datetime'
        :type in_units: str, optional
        :returns:  time point index
        :rtype: {int}
        :raises: NotImplementedError
        """

        if in_units == 'datetime':
            if t < self.dancecard_start_dt:
                raise ValueError('t (%s) is before dancecard start (%s)'%(t.isoformat(),self.dancecard_start_dt.isoformat()))

            return ceil((t - self.dancecard_start_dt).total_seconds() / self.tstep_sec)
        else:
            raise NotImplementedError

    def get_objects_at_ts_pre_tp_indx(self,tp_indx):
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

        if self.mode == 'timepoint':
            raise RuntimeError("this method can't be used in timepoint mode")

        return self.dancecard[tp_indx-1]

    def get_objects_at_ts_indx(self,ts_indx):
        """ get objects stored in the dance card at timestep index
        
        Gets the objects stored in the dance card during the timestep at index ts_indx. 
        :param ts_indx:  index for the timestep of interest
        :type ts_indx: int
        :returns: objects at index
        :rtype: {list}
        """

        if self.mode == 'timepoint':
            raise RuntimeError("this method can't be used in timepoint mode")

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

    @staticmethod
    def get_ts_indx(toi, base_time, tstep_sec):
        '''
        Get the index of the timestep (of duration tstep_sec) in which toi of interest falls

        e.g.
        Dancecard structure:
        |  i=0    |  i=1    |  i=2    |  i=3    ...   |  i=N    |
        | tstep 0 | tstep 1 | tstep 2 | tstep 3 ...   | tstep N |
                       ^
                  *   toi  * (toi anywhere btwn *) --> returns 1
        ^
        base_time

        Uses a floor function, so if toi equals timepoint before timestep, the immediately following timestep index is returned. If toi is almost, but not quite, at the timepoint at the end of a timestep, still returns index of current ts.

        :param datetime toi: time of interest
        :param datetime base_time: first timepoint of dancecard
        :param float tstep_sec: the size of timestep in the dancecard
        :return: the index of timestep containing toi (inclusive of first timepoint, exclusive of trailing tp)
        '''

        toi_sec = (toi - base_time).total_seconds()
        indx = int(floor(toi_sec / tstep_sec))

        return indx

    def get_ts_indx_from_t(self, t, in_units='datetime'):
        """get index of time step containing time t
        
        Get the index of timestep in which input time of interest (t) is contained, using the base time for this dancecard.

        :param t:  the time of interest
        :type t: [specified by in_units]
        :param in_units: type for t, defaults to 'datetime'
        :type in_units: str, optional
        :returns:  time step index
        :rtype: {int}
        :raises: NotImplementedError
        """

        if in_units == 'datetime':
            return Dancecard.get_ts_indx(t, self.dancecard_start_dt, self.tstep_sec)
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
        if self.mode == 'timepoint':
            raise RuntimeError("this method can't be used in timepoint mode")

        if not self.item_type == list:
            raise RuntimeError("this method can't be used if this dancecard stores something other than lists")

        for wind in winds:
            self.add_item_in_interval(wind, wind.start, wind.end)

            # act_start_indx = Dancecard.get_ts_indx(wind.start, self.dancecard_start_dt, self.tstep_sec)
            # act_end_indx = Dancecard.get_ts_indx(wind.end, self.dancecard_start_dt, self.tstep_sec)

            # act_start_indx = max(0, act_start_indx)
            # act_end_indx = min(dancecard_last_indx, act_end_indx)

            # for indx in range(act_start_indx, act_end_indx + 1): # Make sure to include end index
            #     # add object to schedule. Note this is not a deepcopy!

            #     # if list is not initialized yet
            #     if not self.dancecard[indx]:
            #         self.dancecard[indx] = [wind]
            #     else:
            #         self.dancecard[indx].append(wind) 

    def add_item_in_interval(self,item,start,end):

        if self.mode == 'timepoint':
            # round to nearest timepoint in the dancecard. Timestep should not be large enough that this will break things
            # (this is argued from a resource storage perspective) use post tp indx here so that if an activity overlaps at one timepoint with an activity before it, we won't end up overestimating the amount of resource requirement (assumption: two back to back actitivities that overlap at an infinitesimal point do not really overlap)
            start_indx = self.get_tp_indx_post_t(start)
            # use pre tp indx here to be consistent with above.
            end_indx = self.get_tp_indx_pre_t(end)


        elif self.mode == 'timestep':
            dancecard_last_indx = self.num_timesteps - 1

            # todo: this usage of static method is a little stuffy, should update at some point to look like above code for timepoint
            start_indx = Dancecard.get_ts_indx(start, self.dancecard_start_dt, self.tstep_sec)
            end_indx = Dancecard.get_ts_indx(end, self.dancecard_start_dt, self.tstep_sec)

            start_indx = max(0, start_indx)
            end_indx = min(dancecard_last_indx, end_indx)


        for indx in range(start_indx, end_indx + 1): # Make sure to include end index
            # add object to schedule. Note this is not a deepcopy!

            if self.item_type == list:
                # if list is not initialized yet - we have None at this index
                if not self.dancecard[indx]:
                    self.dancecard[indx] = [item]
                else:
                    self.dancecard[indx].append(item) 

            # if it's not a list-based dancecard, just set index equal to item
            else:
                self.dancecard[indx] = item


    def remove_winds_from_dancecard(self, winds, unmodified_yes = False):
        '''
        Remove a set of windows from the dancecard according to their start and stop times. Remove only the objects corresponding to the winds from the dancecard

        :param winds: list of windows to remove. Can be a list with a single element, of course
        :return:
        '''
        if self.mode == 'timepoint':
            raise RuntimeError("this method can't be used in timepoint mode")

        dancecard_last_indx = self.num_timesteps - 1

        for wind in winds:
            act_start_indx = Dancecard.get_ts_indx(wind.start, self.dancecard_start_dt, self.tstep_sec)
            act_end_indx = Dancecard.get_ts_indx(wind.end, self.dancecard_start_dt, self.tstep_sec)

            act_start_indx = max(0, act_start_indx)
            act_end_indx = min(dancecard_last_indx, act_end_indx)

            for indx in range(act_start_indx, act_end_indx + 1): #  make sure to include end index

                try:
                    self.dancecard[indx].remove(wind)  # remove object
                except ValueError:
                    pass  #already been removed