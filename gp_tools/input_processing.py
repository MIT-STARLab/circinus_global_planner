# Contains functionality for turning input data structures into the 
# objects used by the global planner  scheduling module
# 
# 
# @author Kit Kennedy
#

import collections 
from copy import copy, deepcopy

from circinus_tools  import time_tools as tt
from circinus_tools  import  constants as const
from .custom_activity_window import   ObsWindow,  DlnkWindow, XlnkWindow, EclipseWindow
from .schedule_objects  import Dancecard
from .routing_objects import LinkInfo

class GPProcessorIO():
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
        self.num_gs=params['num_gs']

        self.obs_times=params['obs_times']
        self.pl_data_rate=params['pl_data_rate']
        self.targ_ignore_list=params['targ_ignore_list']

        self.dlnk_times=params['dlnk_times']
        self.dlnk_rates=params['dlnk_rates']
        self.min_allowed_dv_dlnk=params['min_allowed_dv_dlnk_Mb']
        self.gs_ignore_list=params['gs_ignore_list']

        self.xlnk_times=params['xlnk_times']
        self.xlnk_rates=params['xlnk_rates']
        self.min_allowed_dv_xlnk=params['min_allowed_dv_xlnk_Mb']

        self.eclipse_times=params['eclipse_times']


    def merge_sat_obs_windows(self,obs_window_list,next_window_uid):
        '''
        Use obs_window_list to create a new list of non-overlapping obs windows, in which each obs activity includes a list of ALL targets being observed at all times.

        :param obs_window_list: the old list of possibly-overlapping obs windows
        :param next_window_uid: unique window ID for unique identification of the windows
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
                    new_obs_window_list.append(ObsWindow(next_window_uid,sat_indx,curr_id_list,sat_target_indx=sat_target_indx,is_urgent=False,start=obs_start,end=obs_end))
                    next_window_uid += 1
                    sat_target_indx += 1

                # refresh target list
                curr_id_list = target_ID_list
                obs_start = copy(obs_end)

        return new_obs_window_list,next_window_uid

    def import_obs_winds( self,next_window_uid=0):
        """  Turn observation times into observation windows

        Parse input data structure to create observation windows. Uses
          indexing  of the input data structure
        
        :returns: [description]
        :rtype: {[type]}
        """
        obs_winds = []
        for sat_indx, all_sat_obs in enumerate(self.obs_times):
            sat_obs_winds = []

            for targ_indx, target_obs in enumerate(all_sat_obs):

                # TODO:  this should operate on ID not on index
                if targ_indx in self.targ_ignore_list:
                    continue

                for obs_indx, obs in enumerate(target_obs):

                    #   convert input date format over to datetime
                    if self.input_date_format == const.MODIFIED_JULIAN_DATE:
                        start =tt.mjd2datetime(obs[0])
                        end=tt.mjd2datetime(obs[1])
                    else:
                        raise NotImplementedError

                    sat_obs_winds.append(ObsWindow(next_window_uid,sat_indx,[targ_indx],sat_target_indx=obs_indx,is_urgent=False,start= start,end= end))
                    next_window_uid+=1

            obs_winds.append(sat_obs_winds)

        for sat_indx in range(self.num_sats):
             new_windows,next_window_uid = self.merge_sat_obs_windows(obs_winds[sat_indx],next_window_uid)
             for wind in new_windows:
                 wind.calc_data_vol(self.pl_data_rate)
             obs_winds[sat_indx] = new_windows

        return obs_winds, next_window_uid

    def import_xlnk_winds( self, next_window_uid=0, sort= True):
        xlink_winds_flat = [[] for i in range(self.num_sats)]
        xlink_winds = [[[] for j in range(self.num_sats)] for i in range(self.num_sats)]
        for sat_indx in range(self.num_sats):
            # xlnk_times  matrix should be symmetrical, so there's no reason to look at lower left  triangle
            for xsat_indx in range(sat_indx+1,self.num_sats):
                xlnk_list = self.xlnk_times[sat_indx][xsat_indx]

                for xlnk_indx, xlnk in enumerate(xlnk_list):

                    #   convert input date format over to datetime
                    if self.input_date_format == const.MODIFIED_JULIAN_DATE:
                        start =tt.mjd2datetime(xlnk[0])
                        end=tt.mjd2datetime(xlnk[1])
                    else:
                        raise NotImplementedError

                    # create a new window 
                    new_wind = XlnkWindow(next_window_uid,sat_indx,xsat_indx,xlnk_indx, start, end)

                    # figure out the data volume for this window
                    xlnk_rates_mat =  self.xlnk_rates[new_wind.sat_indx][new_wind.xsat_indx][new_wind.sat_xsat_indx]
                    new_wind.rates_mat = xlnk_rates_mat
                    new_wind.set_data_vol_and_refresh_times()

                    if new_wind.data_vol >  self.min_allowed_dv_xlnk:
                        #  add to regular matrix
                        xlink_winds[sat_indx][xsat_indx].append(new_wind)
                        # add it to the  flat lists for the sats on both ends of the crosslink. Note that the same object is stored for both, so any modification of the object by one sat modifies it for the other sat as well
                        xlink_winds_flat[sat_indx].append(new_wind)
                        xlink_winds_flat[xsat_indx].append(new_wind)

                    next_window_uid+=1

            # sort the xlink windows for convenience
            if sort:
                xlink_winds_flat[sat_indx].sort(key=lambda x: x.start)


        return  xlink_winds, xlink_winds_flat, next_window_uid

    def import_dlnk_winds( self,next_window_uid=0,sort= True):
        # Import sat dlnk windows
        dlink_winds_flat = []
        dlink_winds = [[[] for j in range(self.num_gs)] for i in range(self.num_sats)]
        for sat_indx, all_sat_dlnk in enumerate( self.dlnk_times):
                        
            sat_dlnk_winds = []

            for gs_indx, dlnk_list in enumerate(all_sat_dlnk):

                # TODO:  this should operate on ID not on index
                if gs_indx in self.gs_ignore_list:
                    continue

                for dlnk_indx, dlnk in enumerate(dlnk_list):

                    #   convert input date format over to datetime
                    if self.input_date_format == const.MODIFIED_JULIAN_DATE:
                        start =tt.mjd2datetime(dlnk[0])
                        end=tt.mjd2datetime(dlnk[1])
                    else:
                        raise NotImplementedError

                    new_wind = DlnkWindow(next_window_uid,sat_indx,gs_indx,dlnk_indx,start, end)

                    new_wind.rates_mat =  self.dlnk_rates[new_wind.sat_indx][new_wind.gs_ID][new_wind.sat_gs_indx]
                    new_wind.set_data_vol_and_refresh_times()

                    if new_wind.data_vol >  self.min_allowed_dv_dlnk:

                        dlink_winds[sat_indx][gs_indx].append (new_wind) 
                        sat_dlnk_winds.append(new_wind)

                    next_window_uid+=1

            # sort the downlink windows for convenience
            if sort:
                sat_dlnk_winds.sort(key=lambda x: x.start)

            dlink_winds_flat.append(sat_dlnk_winds)


        return dlink_winds,dlink_winds_flat, next_window_uid

    def import_eclipse_winds( self,next_window_uid=0):
        """  Turn Eclipse times into eclipse windows

        Parse input data structure to create eclipse windows. Uses
          indexing  of the input data structure
        
        :returns: [description]
        :rtype: {[type]}
        """
        ecl_winds = []
        for sat_indx, ecl_times in enumerate(self.eclipse_times):
            sat_ecl_winds = []

            for ecl_indx, ecl in enumerate(ecl_times):

                #   convert input date format over to datetime
                if self.input_date_format == const.MODIFIED_JULIAN_DATE:
                    start =tt.mjd2datetime(ecl[0])
                    end=tt.mjd2datetime(ecl[1])
                else:
                    raise NotImplementedError

                sat_ecl_winds.append(EclipseWindow(next_window_uid,start= start,end= end))
                next_window_uid+=1

            ecl_winds.append(sat_ecl_winds)

        return ecl_winds, next_window_uid

    def extract_flat_windows( self, routes_flat,  copy_windows= False):
        """ extracts all the activity windows used from a set of routes
        
        Note that if not using copy_windows, the input zwindows objects will be modified
        :param routes_flat:  a flat list of all the routes scheduled
        :type routes_flat: [list(routing_objects.DataRoute)]
        :param copy_windows:  make a copy of the windows within the routes before modifying them, defaults to False
        :type copy_windows: bool, optional
        :returns: three lists, one for observations, one for downlinks and one for cross-links ( each one is a singly nested list with the nesting indexed by sat index),  a dictionary containing link info objects for each dlnk/xlnk window with keys being the dlnk and xlnk windows,  and a similar dictionary containing a list of data route indices for each dlnk/xlnk (specifies which data routes went through which window)
        :rtype: {list(list() by sat_indx),list(list() by sat_indx),list(list() by sat_indx),dict(CommWindow: LinkInfo),dict(CommWindow: list()}
        """

        all_obs = []
        all_xlnk = []
        all_dlnk = []

        #  dictionary of named tuples that contain all of the information used for creating output link info.  keys are the window objects themselves
        link_info_by_wind = {}

        #  dictionary of route indices, using windows as keys
        route_indcs_by_wind  = {}

        def copy_choice(wind):
            if copy_windows:
                return deepcopy(wind)
            else:
                return wind

        for dr_indx, dr in enumerate (routes_flat):
            for wind in dr.route:
                # copy the window before we make any changes to it, if so desired
                wind = copy_choice(wind)

                if type (wind)  == ObsWindow:
                    if wind in all_obs:
                        #  get the matching object ( see note above for all_xlnk)
                        indx = all_obs.index(wind)
                        original_wind = all_obs[indx]
                        #  now we can update the schedule data volume on that original object
                        original_wind.scheduled_data_vol += dr.scheduled_dv
                    else:
                        wind.scheduled_data_vol = dr.scheduled_dv
                        all_obs. append (wind)

                elif type (wind)  == XlnkWindow or type (wind)  == DlnkWindow:

                    if type (wind)  == XlnkWindow:

                        if wind in all_xlnk:
                            #  if we reach here then we know that the hash of wind is in the list, but it's not going to be the same object if we have copied it. do some gymnastics to get the matching object that was already added to the list
                            indx = all_xlnk.index(wind)
                            original_wind = all_xlnk[indx]
                            #  now we can update the schedule data volume on that original object
                            original_wind.scheduled_data_vol += dr.scheduled_dv
                        else:
                            wind.scheduled_data_vol = dr.scheduled_dv
                            all_xlnk. append (wind)

                    if type (wind)  == DlnkWindow:
                        if wind in all_dlnk:
                            #  get the matching object ( see note above for all_xlnk)
                            indx = all_dlnk.index(wind)
                            original_wind = all_dlnk[indx]
                            #  now we can update the schedule data volume on that original object
                            original_wind.scheduled_data_vol += dr.scheduled_dv
                        else:
                            wind.scheduled_data_vol = dr.scheduled_dv
                            all_dlnk. append (wind)

                    #  create link info for this window
                    if not wind in link_info_by_wind.keys (): 
                        link_info_by_wind[wind]  = LinkInfo ([dr.ID], wind.data_vol, dr.scheduled_dv )
                    else:
                        # note: don't have to grab the original window here because we're only using it as an index
                        link_info_by_wind[wind].data_routes.append ( dr.ID) 
                        link_info_by_wind[wind].used_data_vol  +=  dr.scheduled_dv 
                    
                    #  add the route index  for this window,  initializing a list for the window if needed
                    if not wind in route_indcs_by_wind.keys (): 
                        route_indcs_by_wind[wind] = []
                    route_indcs_by_wind[wind].append (dr.ID)


        obs_flat = [[] for k in range( self.num_sats)]
        xlnk_flat = [[] for k in range( self.num_sats)]
        dlnk_flat = [[] for k in range( self.num_sats)]


        for obs in all_obs: 
            obs_flat[obs.sat_indx].append(obs)
        for xlnk in all_xlnk: 
            xlnk_flat[xlnk.sat_indx].append(xlnk)
            xlnk_flat[xlnk.xsat_indx].append(xlnk)
        for dlnk in all_dlnk: 
            dlnk_flat[dlnk.sat_indx].append(dlnk)

        for sat_indx in range ( self.num_sats): 
            obs_flat[sat_indx].sort(key=lambda x: x.start)
            xlnk_flat[sat_indx].sort(key=lambda x: x.start)
            dlnk_flat[sat_indx].sort(key=lambda x: x.start)

        return obs_flat, dlnk_flat, xlnk_flat, link_info_by_wind, route_indcs_by_wind

    def make_sat_history_outputs (self, obs_winds_flat, xlnk_winds_flat, dlnk_winds_flat, link_info_by_wind):
        obs_times_flat = [[] for sat_indx in range ( self.num_sats)]
        obs_locations = [[] for sat_indx in range ( self.num_sats)]
        dlnk_times_flat = [[] for sat_indx in range ( self.num_sats)]
        dlnk_partners = [[] for sat_indx in range ( self.num_sats)]
        dlnk_link_info_history_flat = [[] for sat_indx in range ( self.num_sats)]
        xlnk_times_flat = [[] for sat_indx in range ( self.num_sats)]
        xlnk_partners = [[] for sat_indx in range ( self.num_sats)]
        xlnk_link_info_history_flat = [[] for sat_indx in range ( self.num_sats)]

        for sat_indx in range ( self.num_sats): 
            for wind in obs_winds_flat[sat_indx]:
                #  this  observation window could have multiple targets that it seeing, so we need to separate those out into separate windows
                for target in wind.target_IDs:
                    start_mjd = tt.datetime2mjd ( wind.start)
                    end_mjd = tt.datetime2mjd ( wind.end)
                    obs_times_flat[sat_indx].append ( [start_mjd, end_mjd]) 
                    obs_locations[sat_indx].append ( target)  

        for sat_indx in range ( self.num_sats): 
            for wind in xlnk_winds_flat[sat_indx]:
                # want to filter this so we're not duplicating windows for display in cesium  (though orbit viz filters this too)
                # look at both sat_indx and xsat_indx in case those are not in increasing order
                if wind.sat_indx  > sat_indx or wind.xsat_indx  > sat_indx:
                    start_mjd = tt.datetime2mjd ( wind.start)
                    end_mjd = tt.datetime2mjd ( wind.end)
                    xlnk_times_flat[sat_indx].append ( [start_mjd, end_mjd]) 
                    xlnk_link_info_history_flat[sat_indx].append ( [start_mjd, end_mjd, str (link_info_by_wind[wind])]) 
                    xlnk_partners[sat_indx].append ( wind.xsat_indx)

        for sat_indx in range ( self.num_sats): 
            for wind in dlnk_winds_flat[sat_indx]:
                start_mjd = tt.datetime2mjd ( wind.start)
                end_mjd = tt.datetime2mjd ( wind.end)
                dlnk_times_flat[sat_indx].append ( [start_mjd, end_mjd]) 
                dlnk_link_info_history_flat[sat_indx].append ( [start_mjd, end_mjd, str (link_info_by_wind[wind])]) 
                dlnk_partners[sat_indx].append ( wind.gs_ID)


        # data_history =self.create_data_history( obs_winds_flat,dlnk_winds_flat,xlnk_winds_flat)

        # ordered dictionary so we can preserve order in the output file
        outputs = collections.OrderedDict ()
        outputs['obs_times_flat']  = obs_times_flat 
        outputs['obs_locations']  = obs_locations 
        outputs['dlnk_times_flat']  = dlnk_times_flat 
        outputs['dlnk_partners']  = dlnk_partners
        outputs['dlnk_link_info_history_flat']  = dlnk_link_info_history_flat 
        outputs['xlnk_times_flat']  = xlnk_times_flat 
        outputs['xlnk_partners']  = xlnk_partners 
        outputs['xlnk_link_info_history_flat']  = xlnk_link_info_history_flat 
        # outputs['data_history']  =  data_history
        return outputs


    def create_data_history( self,obs_winds_flat,dlink_winds_flat,xlink_winds_flat):

        #  TODO :   fix this code

        all_wind = [[] for k in range(self.num_sats)]
        for sat_indx in range(self.num_sats):
            for obs in obs_winds_flat[sat_indx]:
                all_wind[sat_indx].append(obs)

            for dlnk in dlink_winds_flat[sat_indx]:
                all_wind[sat_indx].append(dlnk)

            for xlnk in xlink_winds_flat[sat_indx]:
                all_wind[sat_indx].append(xlnk)
                # all_wind[xlnk.xsat_indx].append(xlnk)

            all_wind[sat_indx].sort(key=lambda x: x.start)


        sat_dvs = [0 for i in range(self.num_sats)]
        data_history =  []

        for sat_indx in range(self.num_sats):
            data_history_sat = []

            data_history_sat.append([0,0])

            latest_time_sec = 0

            sat_pkts_state = []
            sat_dv_from_pkts = 0

            for dummy, wind in enumerate(all_wind[sat_indx]):
                old_sat_dv = sat_dvs[sat_indx]
                new_sat_dv = old_sat_dv
                start_time_sec = (wind.start -  self.scenario_start).total_seconds()
                end_time_sec = (wind.end -  self.scenario_start).total_seconds()

                # sanity check
                # add the tstep_sec/2 to give some sub-timestep wiggle room - sometimes the times calculated for the windows aren't quite on the timestep border
                print ( start_time_sec)
                print ( all_wind[sat_indx])

                if start_time_sec < (latest_time_sec - self.tstep_sec/2):
                    raise Exception ('create_data_history: problem in data history construction, times wrong')

                if wind.data_vol > wind.unmodified_data_vol or wind.remaining_data_vol < 0:
                    raise Exception ('create_data_history: problem in data history construction, data vol wrong')

                if type(wind) == ObsWindow:
                    dv_created = wind.collected_data_vol
                    new_sat_dv +=  dv_created

                    if dv_created != sum(pkt.data_vol for pkt in wind.data_pkts):
                        raise Exception ('create_data_history: problem 3')

                    for pkt in wind.data_pkts:
                        sat_pkts_state.append(pkt)

                elif type(wind) == DlnkWindow:
                    dv_dlnked = wind.routed_data_vol
                    new_sat_dv -= dv_dlnked

                    if dv_dlnked != sum(pkt.data_vol for pkt in wind.data_pkts):
                        raise Exception ('create_data_history: problem 4')

                    for pkt in wind.data_pkts:
                        sat_pkts_state.remove(pkt)

                elif type(wind) == XlnkWindow:
                    for sat_indx_to in wind.routed_data_vol_to_sat_indx.keys():
                        dv_xlnked = wind.routed_data_vol_to_sat_indx[sat_indx_to]
                        pkts_xlnked = wind.routed_pkts_to_sat_indx[sat_indx_to]

                        if dv_xlnked != sum(pkt.data_vol for pkt in pkts_xlnked):
                            raise Exception ('create_data_history: problem 5')

                        if sat_indx_to == sat_indx:
                            new_sat_dv +=  dv_xlnked
                            for pkt in pkts_xlnked:
                                sat_pkts_state.append(pkt)
                        else:
                            new_sat_dv -= dv_xlnked
                            for pkt in pkts_xlnked:
                                sat_pkts_state.remove(pkt)

                sat_dv_from_pkts = sum(pkt.data_vol for pkt in sat_pkts_state)

                if new_sat_dv < -0.9:
                    print('create_data_history: sat data volume < 0')
                    print('new_sat_dv: '+str(new_sat_dv))
                    print('sat_dv_from_pkts: '+str(sat_dv_from_pkts))
                    # 1/0

                data_history_sat.append([start_time_sec, old_sat_dv])
                data_history_sat.append([end_time_sec, new_sat_dv])

                sat_dvs[sat_indx] = new_sat_dv
                latest_time_sec = max(end_time_sec,latest_time_sec)  # use this max in case for some reason the end times of the windows are bouncing around.

            data_history.append(data_history_sat)

            if len(sat_pkts_state) > 0:
                raise Exception ('create_data_history: should not have any packets left in here')

        return data_history