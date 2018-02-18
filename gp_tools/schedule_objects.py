from copy import copy, deepcopy
from datetime import timedelta

class Dancecard(object):
    def __init__(self, dancecard_start, dancecard_end, tstep_sec):
        '''
        An entity schedule object maintains a schedule data structure that keeps track of the activities that an entity (e.g. satellite, ground station) has committed to perform, and the start and stop times for those acts

        Dancecard structure:
        |  i=0    |  i=1    |  i=2    |  i=3    ...   |  i=N    |
        | tstep 0 | tstep 1 | tstep 2 | tstep 3 ...   | tstep N |
        ^                             ^
        base_time                *    t    *

        :param scheduling_start: the start of the scheduling window (generally same as scenario)
        :param scheduling_end: the end of the scheduling window (generally same as scenario)
        :param tstep_sec: time step in seconds in overall scenario
        '''

        total_duration = (dancecard_end - dancecard_start).total_seconds()
        # each index in dancecard represents a chunk of time, e.g. [1] is from start_time+tstep_sec*1 to start_time+tstep_sec*2
        num_timesteps = int(total_duration / tstep_sec)

        # this "dancecard" stores the objects for the activities that are being executed at each timepoint (each index). Only one activity is possible at a given time
        self.dancecard = [[] for i in range(num_timesteps)]

        self.dancecard_start = dancecard_start
        self.dancecard_end = dancecard_end
        self.tstep_sec = tstep_sec
        self.tstep_td = timedelta(seconds=tstep_sec)
        self.num_timesteps = num_timesteps

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

# class EntitySchedule(Dancecard):
#     def __init__(self, scheduling_start, scheduling_end, tstep_sec, max_num_concurrent_acts):
#         '''
#         An entity schedule object maintains a schedule data structure that keeps track of the activities that an entity (e.g. satellite, ground station) has committed to perform, storing those activities as objects within a dancecard. Start and end times of activities can found by counting indices in the dancecard

#         :param scheduling_start: the start of the scheduling window (generally same as scenario)
#         :param scheduling_end: the end of the scheduling window (generally same as scenario)
#         :param tstep_sec: time step in seconds in overall scenario
#         '''

#         self.max_num_concurrent_acts = max_num_concurrent_acts

#         self._indx_overlap_count = 0

#         super(EntitySchedule, self).__init__(scheduling_start, scheduling_end, tstep_sec)

#     def addWindsToSchedule(self, winds):
#         '''
#         Add a set of windows to the sat's schedule according to their start and stop times, adding to whatever activities are already scheduled, if there's still availability.

#         :param winds: list of windows to add. Can be a list with a single element, of course
#         :return: boolean that is True if any of the winds was successfully added to the schedule, even only at one point in time
#         '''
#         added_to_schedule = False

#         dancecard_last_indx = len(self.dancecard)-1

#         for wind in winds:
#             act_start_indx = Dancecard.getPostIndexDancecard(wind.start, self.dancecard_start, self.tstep_sec)
#             act_end_indx = Dancecard.getPostIndexDancecard(wind.end, self.dancecard_start, self.tstep_sec) - 1 # -1 because last time slice is before wind.end

#             act_start_indx = max(0, act_start_indx)
#             act_end_indx = min(dancecard_last_indx,act_end_indx)

#             for indx in xrange(act_start_indx, act_end_indx + 1): # go all the way up to end indx
#                 if wind in self.dancecard[indx]:
#                     pass  # the object has already been scheduled at this instant in time
#                 elif len(self.dancecard[indx]) == self.max_num_concurrent_acts:
#                     print 'addWindsToSchedule: oh god, schedule overlap!'
#                     self._indx_overlap_count += 1
#                 else:
#                     self.dancecard[indx].append(wind)  # add object to schedule. Note this is not a deepcopy!
#                     added_to_schedule = True

#         return added_to_schedule


#     def extractDeconflictedWindows(self, winds, window_id):
#         '''
#         Takes a list of windows and figures out when these windows DO NOT overlap any already-scheduled activities for this SatSchedule. Returns a list of new windows with no overlap. Note that the objects in the returned list have been deepcopied for safety and simplicity of implementation. Because of deepcopy, this function is potentially expensive to run. Also the returned windows are not filtered for length (should be >= 1 timestep) nor are datavolumes updated.

#         :param ActivityWindow[N] winds: list of windows to be deconflicted. These windows are assumed to be non-overlapping
#         :param window_id: window ID value to start at for assigning ID to new copies of windows
#         :return: ActivityWindow[N+M] deconflict_winds: new list of windows created from winds, but having no overlap with any currently scheduled activities for this satellite (M is number of new windows created from splits)
#         '''

#         # make this a bunch of Nones, not a bunch of empty lists, because of assumptions winds are non-overlapping
#         temp_dancecard = [None for i in xrange(self.num_timesteps)]

#         dancecard_last_indx = len(self.dancecard) - 1

#         #######
#         # Add each wind to the temp dancecard, NOT overwriting what's already there
#         #######

#         for wind in winds:
#             act_start_indx = Dancecard.getPostIndexDancecard(wind.start, self.dancecard_start, self.tstep_sec)
#             act_end_indx = Dancecard.getPostIndexDancecard(wind.end, self.dancecard_start,self.tstep_sec) - 1  # -1 because last time slice is before wind.end

#             act_start_indx = max(0, act_start_indx)
#             act_end_indx = min(dancecard_last_indx, act_end_indx)

#             for indx in xrange(act_start_indx, act_end_indx + 1): # go all the way up to end indx
#                 if len(self.dancecard[indx]) < self.max_num_concurrent_acts:  # if still scheduling availability
#                     temp_dancecard[indx] = wind  # add object to schedule. Note this is not a deepcopy!

#         #######
#         # Add each wind to the temp dancecard, NOT overwriting what's already there
#         #######

#         # note: temp_dancecard only contains winds objects at non-conflicting times. Need to go back through and recreate windows, because originals might have been split

#         last_wind_obj = temp_dancecard[0]
#         wind_start = copy(self.dancecard_start)
#         deconflict_winds = []
#         for indx, wind_obj in enumerate(temp_dancecard):

#             # check if object has changed
#             if not wind_obj == last_wind_obj:
#                 wind_end = self.dancecard_start + timedelta(seconds=indx * self.tstep_sec)

#                 if last_wind_obj:  # if it's not None
#                     # need to update timing information for window - that could have changed due to conflicts with scheduled activities.
#                     # deepcopy the current object. Have to do this because we don't know if this object has actually been split into multiple intervals
#                     wind_obj_to_add = deepcopy(last_wind_obj)
#                     wind_obj_to_add.start = wind_start
#                     wind_obj_to_add.end = wind_end
#                     wind_obj_to_add.refreshDuration()

#                     # window_ID attribute is useful for human-recognizing different windows. Update if it's there.
#                     if hasattr(wind_obj_to_add, 'window_ID'):
#                         wind_obj_to_add.window_ID = window_id
#                         window_id += 1

#                     deconflict_winds.append(wind_obj_to_add)

#                 # refresh target list
#                 last_wind_obj = wind_obj
#                 wind_start = copy(wind_end)

#         return deconflict_winds, window_id

#     def isAvailable(self,t_indx):
#         '''
#         See if this moment in time in the schedule is open for another activity

#         :param t_indx: index into the schedule dancecard
#         :return:
#         '''

#         time_list = self.dancecard[t_indx]
#         if len(time_list) >= self.max_num_concurrent_acts:
#             return False

#         return True

#     def scheduleMatches(self,t_indx,act):
#         '''
#         See if act is already scheduled at this moment in time

#         :param t_indx: index into the schedule dancecard
#         :param act:
#         :return:
#         '''

#         time_list = self.dancecard[t_indx]
#         if act in time_list:
#             return True

#         return False

# class SatSchedule(EntitySchedule):
#     def __init__(self, sat_indx, scheduling_start, scheduling_end, tstep_sec):
#         '''
#         A sat schedule object maintains a schedule data structure that keeps track of the activities that a satellite has committed to perform, and the start and stop times for those acts

#         :param sat_indx: index of satellite owning this schedule
#         :param scheduling_start: the start of the scheduling window (generally same as scenario)
#         :param scheduling_end: the end of the scheduling window (generally same as scenario)
#         :param tstep_sec: time step in seconds in overall scenario
#         '''

#         self.sat_indx = sat_indx

#         num_concurrent_acts = 1

#         super(SatSchedule, self).__init__(scheduling_start, scheduling_end, tstep_sec, num_concurrent_acts)

#     # @staticmethod
#     # def removeFullScheduleFromAvailabilityDancecard(schedule,avail_dancecards,num_sats):
#     #     '''
#     #
#     #
#     #     :param SatSchedule schedule: the sat schedule that we're using to
#     #     :param Dancecard[num_sats] avail_dancecards:
#     #     :return:
#     #     '''
#     #
#     #     sat_avail_dancecard = avail_dancecards[schedule.sat_indx]
#     #
#     #     # check that schedule and dancecard inputs actually cover same time period
#     #     if not (schedule.dancecard_start == sat_avail_dancecard.dancecard_start and  schedule.dancecard_end == sat_avail_dancecard.dancecard_end):
#     #         print 'schedule and dancecard must have same timesteps!'
#     #         1/0
#     #
#     #
#     #     for indx, time_list in enumerate(schedule.dancecard):
#     #
#     #         # if this particular index in the schedule is already full with activities, remove from availability dancecard
#     #         if len(time_list)==schedule.max_num_concurrent_acts:
#     #             sat_avail_dancecard[indx] = []

#                 # overlap_acts_to_remove = sat_avail_dancecard.dancecard[indx]

#             # # it's unfortunate, but we have to take all overlapping activities we found during this time and search for them across the dancecards and remove them. This comes about because overlap_acts can be xlnks. There are two ends on the xlnk - if we just remove the one end at act.sat_indx (and act.xsat_indx) then we'll have a loose end at another sat index where that sat thinks it can still crosslink with this sat
#             # for overlap_act in overlap_acts_to_remove:
#             #     for other_sat_indx in xrange(num_sats):
#             #
#             #         other_sat_dancecard = avail_dancecards[other_sat_indx].dancecard[indx]
#             #         if overlap_act in other_sat_dancecard:
#             #             other_sat_dancecard.remove(overlap_act) # todo: is using remove here super inefficient?


# class GSSchedule(EntitySchedule):
#     def __init__(self, gs_indx, gs_name, scheduling_start, scheduling_end, tstep_sec):
#         '''
#         A ground station schedule object maintains a schedule data structure that keeps track of the satellites with which a ground station is currently engaged in a downlink, and the start and stop times for those acts

#         :param gs_indx: index of ground station owning this schedule
#         :param gs_name: name of ground station owning this schedule
#         :param scheduling_start: the start of the scheduling window (generally same as scenario)
#         :param scheduling_end: the end of the scheduling window (generally same as scenario)
#         :param tstep_sec: time step in seconds in overall scenario
#         '''

#         self.gs_indx = gs_indx
#         self.gs_name = gs_name
#         num_concurrent_acts = 1

#         super(GSSchedule, self).__init__(scheduling_start, scheduling_end, tstep_sec, num_concurrent_acts)

