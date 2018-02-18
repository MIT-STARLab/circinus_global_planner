# Contains functionality for turning input data structures into the 
# objects used by the global planner  scheduling module
# 
# 
# @author Kit Kennedy
#

import collections 
from copy import copy

from circinus_tools  import time_tools as tt
from circinus_tools  import  constants as const
from .custom_activity_window import   ObsWindow,  DlnkWindow, XlnkWindow
from .schedule_objects  import Dancecard

class GPInputProcessor():
    """docstring for GPInputProcessor"""
    def __init__(self,params):
        #  assume modified Julian date for now. todo: make this a parameter
        self.input_date_format = const.MODIFIED_JULIAN_DATE

        #   note: we're implicitly testing for the presence of these parameters  below

        #  common parameters used for processing  multiple types of input windows
        # scenario_start: datetime storing the start of the overall simulation period
        # scenario_end: datetime storing the end of the overall simulation period
        # tstep_sec: timestep for sim, in seconds
        self.scenario_start=params['start_utc_dt']
        self.scenario_end=params['end_utc_dt']
        self.tstep_sec=params['timestep_s']
        self.num_sats=params['num_sats']

        self.obs_times=params['obs_times']
        self.pl_data_rate=params['pl_data_rate']
        self.targ_ignore_list=params['targ_ignore_list']

        self.dlnk_times=params['dlnk_times']
        self.dlnk_rates=params['dlnk_rates']
        self.min_allowed_dv_dlnk=params['min_allowed_dv_dlnk']
        self.gs_ignore_list=params['gs_ignore_list']

        self.xlnk_times=params['xlnk_times']
        self.xlnk_rates=params['xlnk_rates']
        self.min_allowed_dv_xlnk=params['min_allowed_dv_xlnk']



    def merge_sat_obs_windows(self,obs_window_list,window_id):
        '''
        Use obs_window_list to create a new list of non-overlapping obs windows, in which each obs activity includes a list of ALL targets being observed at all times.

        :param obs_window_list: the old list of possibly-overlapping obs windows
        :param window_id: latest window ID for obs windows, merely for convenience of IDing the windows
        :return: new obs list with non-overlapping obs events that replaces input obs_window_list
        '''

        dc_target_IDs = Dancecard(self.scenario_start, self.scenario_end, self.tstep_sec)

        sat_indx = obs_window_list[0].sat_indx

        # populate the dancecard with the obs targets at each time
        for obs in obs_window_list:
            obs_start = obs.start
            obs_end = obs.end
            if obs_start > self.scenario_end:
                continue
            elif obs_end > self.scenario_end:
                obs_end = self.scenario_end

            obs_start_indx = dc_target_IDs.get_post_index_by_time(obs_start)
            obs_end_indx = dc_target_IDs.get_post_index_by_time(obs_end) - 1 # -1 because last time slice is before wind.end

            for indx in range(obs_start_indx,obs_end_indx+1):
                dc_target_IDs.dancecard[indx] += obs.target_IDs

        curr_id_list = dc_target_IDs.dancecard[0]
        obs_start = copy(self.scenario_start)
        new_obs_window_list = []
        sat_target_indx = 0
        for indx, target_ID_list in enumerate(dc_target_IDs.dancecard):

            # check if lists have the all the same elements in each other (order not important)
            if not collections.Counter(target_ID_list) == collections.Counter(curr_id_list):
                # obs_end = self.scenario_start + timedelta(seconds=(indx+1)*self.tstep_sec)
                obs_end = dc_target_IDs.get_pre_time_by_index(indx+1)  # plus one on index because we're looking for abs time after last timestep

                # todo: it's should probably add filtering for minimum length observations here

                if len(curr_id_list) > 0:  # if it's not empty
                    # create a new observation based on the previous set of targets
                    new_obs_window_list.append(ObsWindow(window_id,sat_indx,curr_id_list,sat_target_indx=sat_target_indx,is_urgent=False,start=obs_start,end=obs_end))
                    window_id += 1
                    sat_target_indx += 1

                # refresh target list
                curr_id_list = target_ID_list
                obs_start = copy(obs_end)

        return new_obs_window_list

    def import_obs_winds( self):
        """  Turn observation times into observation windows

        Parse input data structure to create observation windows. Uses
          indexing  of the input data structure
        
        :returns: [description]
        :rtype: {[type]}
        """
        obs_winds = []
        obs_window_id = 0
        for sat_indx, all_sat_obs in enumerate(self.obs_times):
            sat_obs_winds = []

            for targ_indx, target_obs in enumerate(all_sat_obs):

                if targ_indx in self.targ_ignore_list:
                    continue

                for obs_indx, obs in enumerate(target_obs):

                    #   convert input date format over to datetime
                    if self.input_date_format == const.MODIFIED_JULIAN_DATE:
                        start =tt.mjd2datetime(obs[0])
                        end=tt.mjd2datetime(obs[1])
                    else:
                        raise NotImplementedError

                    sat_obs_winds.append(ObsWindow(obs_window_id,sat_indx,[targ_indx],sat_target_indx=obs_indx,is_urgent=False,start= start,end= end))
                    obs_window_id+=1

            obs_winds.append(sat_obs_winds)

        for sat_indx in range(self.num_sats):
             new_windows = self.merge_sat_obs_windows(obs_winds[sat_indx],obs_window_id)
             for wind in new_windows:
                 wind.calc_data_vol(self.pl_data_rate)
             obs_winds[sat_indx] = new_windows

        return obs_winds, obs_window_id

    def import_xlnk_winds( self, sort= True):
        xlink_winds = [[] for j in range(self.num_sats)]
        xlnk_window_id = 0
        for sat_indx in range(self.num_sats):
            for xsat_indx in range(sat_indx+1,self.num_sats):
                xlnk_list = self.xlnk_times[sat_indx][xsat_indx]

                for xlnk_indx, xlnk in enumerate(xlnk_list):

                    #   convert input date format over to datetime
                    if self.input_date_format == const.MODIFIED_JULIAN_DATE:
                        start =tt.mjd2datetime(xlnk[0])
                        end=tt.mjd2datetime(xlnk[1])
                    else:
                        raise NotImplementedError

                    # create a new window and add it to the lists for the sats on both ends of the crosslink. Note that the same object is stored for both, so any modification of the object by one sat modifies it for the other sat as well
                    new_wind = XlnkWindow(xlnk_window_id,sat_indx,xsat_indx,xlnk_indx, start, end)
                    xlink_winds[sat_indx].append(new_wind)
                    xlink_winds[xsat_indx].append(new_wind)
                    xlnk_window_id+=1

            # sort the xlink windows for convenience
            if sort:
                xlink_winds[sat_indx].sort(key=lambda x: x.start)

        # calculate data volumes for xlnk windows
        windows_to_keep = []
        for window_list in xlink_winds:
            windows_to_keep_temp = []
            for wind in window_list:

                # half of these xlnk windows are duplicated due to symmetry, so check to see if data volume hasn't been calced yet before trying
                if wind.data_vol == const.UNASSIGNED:
                    xlnk_rates_mat =  self.xlnk_rates[wind.sat_indx][wind.xsat_indx][wind.sat_xsat_indx]

                    wind.rates_mat = xlnk_rates_mat
                    wind.set_data_vol_and_refresh_times()

                    if wind.data_vol >  self.min_allowed_dv_xlnk:
                        windows_to_keep_temp.append(wind)

            windows_to_keep.append(windows_to_keep_temp)

        xlink_winds = windows_to_keep

        return  xlink_winds, xlnk_window_id

    def import_dlnk_winds( self,sort= True):
        # Import sat dlnk windows
        dlink_winds = []
        dlnk_window_id = 0
        for sat_indx, all_sat_dlnk in enumerate( self.dlnk_times):
            sat_dlnk_winds = []

            for gs_indx, dlnk_list in enumerate(all_sat_dlnk):

                if gs_indx in self.gs_ignore_list:
                    continue

                for dlnk_indx, dlnk in enumerate(dlnk_list):

                    #   convert input date format over to datetime
                    if self.input_date_format == const.MODIFIED_JULIAN_DATE:
                        start =tt.mjd2datetime(dlnk[0])
                        end=tt.mjd2datetime(dlnk[1])
                    else:
                        raise NotImplementedError

                    new_wind = DlnkWindow(dlnk_window_id,sat_indx,gs_indx,dlnk_indx,start, end)
                    sat_dlnk_winds.append(new_wind)
                    dlnk_window_id+=1

            # sort the downlink windows for convenience
            if sort:
                sat_dlnk_winds.sort(key=lambda x: x.start)

            dlink_winds.append(sat_dlnk_winds)

        # calculate data volumes for dlnk windows
        windows_to_keep = []
        for window_list in dlink_winds:
            windows_to_keep_temp = []
            for wind in window_list:

                wind.rates_mat =  self.dlnk_rates[wind.sat_indx][wind.gs_ID][wind.sat_gs_indx]
                wind.set_data_vol_and_refresh_times()

                if wind.data_vol >  self.min_allowed_dv_dlnk:
                    windows_to_keep_temp.append(wind)

            windows_to_keep.append(windows_to_keep_temp)

        dlink_winds = windows_to_keep

        return dlink_winds, dlnk_window_id