#  contains code for assessing metrics of global planner output
#
# @author Kit Kennedy
#
#  note that a path is the same as a route.

import time
from datetime import datetime, timedelta
import numpy as np

from circinus_tools  import time_tools as tt
from .custom_activity_window import   ObsWindow,  DlnkWindow, XlnkWindow
from .routing_objects import DataRoute, DataMultiRoute

class GPMetrics():
    """docstring for GPMetrics"""

    ALLOWED_DR_TYPES = [DataRoute, DataMultiRoute]

    def __init__(self, gp_params):
        """initializes based on parameters
        
        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """
        scenario_params = gp_params['gp_orbit_prop_params']['scenario_params']
        sat_params = gp_params['gp_orbit_prop_params']['sat_params']
        gp_inst_params = gp_params['gp_instance_params']['metrics_params']
        obs_params = gp_params['gp_orbit_prop_params']['obs_params']
        gp_general_other_params = gp_params['gp_general_params']['other_params']
        metrics_params = gp_params['gp_general_params']['metrics_params']
        plot_params = gp_params['gp_general_params']['plot_params']
        as_params = gp_params['gp_general_params']['activity_scheduling_params']

        self.latency_params = gp_params['gp_general_params']['other_params']['latency_calculation']
        self.scenario_start_dt  = scenario_params['start_utc_dt']
        self.met_start_dt  = tt.iso_string_to_dt (gp_inst_params['start_utc'])
        self.met_end_dt  = tt.iso_string_to_dt (gp_inst_params['end_utc'])
        self.num_sats=sat_params['num_sats']
        self.num_targ = obs_params['num_targets']
        self.all_targ_IDs = [targ['id'] for targ in obs_params['targets']]
        self.targ_id_ignore_list = gp_general_other_params['targ_id_ignore_list']

        self.min_as_route_dv = as_params['min_as_route_dv_Mb']
        self.aoi_units = metrics_params['aoi_units']
        self.overlap_count_option = metrics_params['overlap_count_option']
        self.window_overlap_option = metrics_params['window_overlap_option']
        self.overlap_window_td = timedelta(metrics_params['overlap_window_s'])
        self.overlap_window_s = metrics_params['overlap_window_s']
        self.overlap_include_dlnks = metrics_params['overlap_include_dlnks']
        self.overlap_remove_same_obs = metrics_params['overlap_remove_same_obs']
        self.aoi_plot_t_units=plot_params['time_units']

        self.power_params = sat_params['power_params_sorted']
        self.sats_emin_Wh = [p_params['battery_storage_Wh']['e_min'][p_params['battery_option']] for p_params in self.power_params]
        self.sats_emax_Wh = [p_params['battery_storage_Wh']['e_max'][p_params['battery_option']] for p_params in self.power_params]

        # the amount by which the minimum data volume is allowed to be lower than self.min_as_route_dv
        self.min_as_obs_dv_slop = self.min_as_route_dv*0.01

        # if two downlink times are within this number of seconds, then they are counted as being at the same time for the purposes of AoI calculation
        self.dlnk_same_time_slop_s = scenario_params['timestep_s'] - 1

    def assess_dv_all_routes(self, routes,  verbose  = False):
        stats = {}

        dvs = []
        for dr in routes:
            dvs.append(dr.scheduled_dv)

        valid = len(dvs) > 0

        stats['total_dv'] = sum(dvs) if valid else 0
        stats['ave_rt_dv'] = np.mean(dvs) if valid else 0
        stats['min_rt_dv'] = np.min(dvs) if valid else 0
        stats['max_rt_dv'] = np.max(dvs) if valid else 0

        if verbose:
            print('data volume for routes')
            print("%s: %f"%('total_dv',stats['total_dv']))
            print("%s: %f"%('ave_rt_dv',stats['ave_rt_dv']))
            print("%s: %f"%('min_rt_dv',stats['min_rt_dv']))
            print("%s: %f"%('max_rt_dv',stats['max_rt_dv']))

        return stats

    @staticmethod
    def get_routes_by_obs(routes):
        rts_by_obs = {}
        for dr in routes:
            obs = dr.get_obs()
            if not obs in rts_by_obs.keys ():
                rts_by_obs[obs] = []
                
            rts_by_obs[obs].append (dr)

        return rts_by_obs

    def assess_dv_by_obs(self, rs_routes_by_obs, sched_routes, verbose=False):
        """ assess data volume down linked as grouped by observation
        
        pretty straightforward, we figure out all the data routes for each observation, and calculate statistics on data volumes for each observation
        :param routes: [description]
        :type routes: [type]
        :param verbose: [description], defaults to False
        :type verbose: bool, optional
        :returns: [description]
        :rtype: {[type]}
        """

        stats = {}

        sched_rts_by_obs = self.get_routes_by_obs (sched_routes)

        rs_dvs_by_obs =  {}
        sched_dvs_by_obs =  {}
        num_sched_obs = 0
        num_rs_obs_dv_not_zero = 0
        for obs in rs_routes_by_obs.keys ():
            rs_dvs_by_obs[obs] = min(obs.data_vol,sum (rt.data_vol for rt in rs_routes_by_obs[obs]))
            if rs_dvs_by_obs[obs] > 0:
                num_rs_obs_dv_not_zero += 1
            if obs in sched_rts_by_obs.keys():
                sched_dvs_by_obs[obs] = sum (rt.scheduled_dv for rt in sched_rts_by_obs[obs])
                num_sched_obs +=1
            else:
                sched_dvs_by_obs[obs] = 0


        rs_dvs = [dv for dv in rs_dvs_by_obs. values ()]
        sched_dvs = [dv for dv in sched_dvs_by_obs. values ()]
        
        valid = len(rs_dvs) > 0

        stats['num_obs_rs'] = len(rs_dvs)
        stats['num_obs_rs_pos_dv'] = num_rs_obs_dv_not_zero
        stats['num_obs_sched'] = num_sched_obs
        stats['total_collectible_dv'] = sum(rs_dvs) if valid else 0
        stats['total_sched_dv'] = sum(sched_dvs) if valid else 0
        stats['ave_obs_dv_rs'] = np.mean(rs_dvs) if valid else 0
        stats['ave_obs_dv_sched'] = np.mean(sched_dvs) if valid else 0
        stats['std_obs_dv_rs'] = np.std(rs_dvs) if valid else 0
        stats['std_obs_dv_sched'] = np.std(sched_dvs) if valid else 0
        stats['min_obs_dv_rs'] = np.min(rs_dvs) if valid else 0
        stats['min_obs_dv_sched'] = np.min(sched_dvs) if valid else 0
        stats['max_obs_dv_rs'] = np.max(rs_dvs) if valid else 0
        stats['max_obs_dv_sched'] = np.max(sched_dvs) if valid else 0

        stats['rs_dvs_by_obs'] = rs_dvs_by_obs

        if verbose:
            print('------------------------------')
            print('data volume by observation')
            print("%s: %f"%('num_obs_rs',stats['num_obs_rs']))
            print("%s: \t\t\t %f"%('num_obs_rs_pos_dv',stats['num_obs_rs_pos_dv']))
            print("%s: \t\t\t\t %f"%('num_obs_sched',stats['num_obs_sched']))
            print("%s: \t\t\t %f"%('total_collectible_dv',stats['total_collectible_dv']))
            print("%s: \t\t\t %f"%('total_sched_dv',stats['total_sched_dv']))
            print("%s: %f"%('ave_obs_dv_rs',stats['ave_obs_dv_rs']))
            print("%s: %f"%('std_obs_dv_rs',stats['std_obs_dv_rs']))
            print("%s: %f"%('min_obs_dv_rs',stats['min_obs_dv_rs']))
            print("%s: %f"%('max_obs_dv_rs',stats['max_obs_dv_rs']))
            print("%s: %f"%('ave_obs_dv_sched',stats['ave_obs_dv_sched']))
            print("%s: %f"%('std_obs_dv_sched',stats['std_obs_dv_sched']))
            print("%s: %f"%('min_obs_dv_sched',stats['min_obs_dv_sched']))
            print("%s: %f"%('max_obs_dv_sched',stats['max_obs_dv_sched']))

            # for obs, dv in dvs_by_obs.items ():
            #     print("%s: %f"%(obs,dv))

        return stats

    def assess_latency_all_routes(self, routes,verbose  = False):
        """ assess latency for all routes
        
        pretty straightforward, just calculate latency for every route and then do statistics on that
        :param routes: [description]
        :type routes: [type]
        :param verbose: [description], defaults to False
        :type verbose: bool, optional
        :returns: [description]
        :rtype: {[type]}
        """
        stats = {}

        latencies = []
        for dr in routes:
            latencies.append(
                dr.get_latency(
                    'minutes',
                    obs_option = self.latency_params['obs'], 
                    dlnk_option = self.latency_params['dlnk']
                )
            )

        valid = len(latencies) > 0

        stats['ave_lat_mins'] = np.mean(latencies) if valid else None
        stats['min_lat_mins'] = np.min(latencies) if valid else None
        stats['max_lat_mins'] = np.max(latencies) if valid else None

        if verbose and valid:
            print('------------------------------')
            print('latency for routes')
            print("%s: %f"%('ave_lat_mins',stats['ave_lat_mins']))
            print("%s: %f"%('min_lat_mins',stats['min_lat_mins']))
            print("%s: %f"%('max_lat_mins',stats['max_lat_mins']))


        return stats

    def assess_latency_by_obs(self, rs_routes_by_obs,sched_routes, verbose=False):
        """ assess downlink latency as grouped by observation
        
        less straightforward than latency by route. First we group by observation,  then we find out how long it took to downlink the first minimum desired amount of data for each observation. based on how long this took we determin the latency of downlink for the observation.
        :param sched_routes: [description]
        :type sched_routes: [type]
        :param verbose: [description], defaults to False
        :type verbose: bool, optional
        :returns: [description]
        :rtype: {[type]}
        """

        stats = {}

        #  for selected routes
        intial_lat_by_obs_rs =  {}
        for obs, rts in rs_routes_by_obs.items ():
            # start, center, end...whichever we're using for the latency calculation
            time_option = self.latency_params['dlnk']

            #  want to sort these by earliest time so that we favor earlier downlinks
            rts.sort (key=lambda rt: getattr(rt.get_dlnk(),time_option))

            #  figure out the latency for the initial minimum DV downlink
            #  have to accumulate data volume because route selection minimum data volume might be less than that for activity scheduling
            cum_dv = 0
            for dr in rts:
                cum_dv += dr.data_vol
                
                #  if we have reached our minimum required data volume amount to deem the observation downlinked for the purposes of latency calculation...
                if cum_dv >= self.min_as_route_dv - self.min_as_obs_dv_slop :

                    intial_lat_by_obs_rs[obs] = dr.get_latency(
                        'minutes',
                        obs_option = self.latency_params['obs'], 
                        dlnk_option = self.latency_params['dlnk']
                    )

                    #  break so that we don't continue considering the rest of the data volume
                    break

        #  for scheduled routes
        sched_rts_by_obs = self.get_routes_by_obs (sched_routes)
        intial_lat_by_obs =  {}
        final_lat_by_obs =  {}
        for obs, rts in sched_rts_by_obs.items ():
            # start, center, end...whichever we're using for the latency calculation
            time_option = self.latency_params['dlnk']

            #  want to sort these by earliest time so that we favor earlier downlinks
            rts.sort (key=lambda rt: getattr(rt.get_dlnk(),time_option))

            #  figure out the latency for the first route that got downlinked.
            # sanity check that its scheduled data volume meets the minimum requirement
            assert(rts[0].scheduled_dv >= self.min_as_route_dv - self.min_as_obs_dv_slop)
            intial_lat_by_obs[obs] = rts[0].get_latency(
                'minutes',
                obs_option = self.latency_params['obs'], 
                dlnk_option = self.latency_params['dlnk']
            )

            # figure out the latency for downlink of all observation data that we chose to downlink
            final_lat_by_obs[obs] = rts[-1].get_latency(
                'minutes',
                obs_option = self.latency_params['obs'], 
                dlnk_option = self.latency_params['dlnk']
            )

        i_lats_rs = [lat for lat in intial_lat_by_obs_rs. values ()]
        i_lats = [lat for lat in intial_lat_by_obs. values ()]
        f_lats = [lat for lat in final_lat_by_obs. values ()]
        
        i_valid = len(i_lats) > 0
        f_valid = len(f_lats) > 0

        # from circinus_tools import debug_tools
        # debug_tools.debug_breakpt()

        #  note that if center times are not used  with both the observation and downlink for calculating latency, then the route selection and scheduled the numbers might differ. this is because the scheduled numbers reflect the updated start and end time for the Windows
        stats['ave_obs_initial_lat_rs'] = np.mean(i_lats_rs) if i_valid else 0
        stats['std_obs_initial_lat_rs'] = np.std(i_lats_rs) if i_valid else 0
        stats['min_obs_initial_lat_rs'] = np.min(i_lats_rs) if i_valid else 0
        stats['max_obs_initial_lat_rs'] = np.max(i_lats_rs) if i_valid else 0
        stats['ave_obs_initial_lat'] = np.mean(i_lats) if i_valid else 0
        stats['std_obs_initial_lat'] = np.std(i_lats) if i_valid else 0
        stats['min_obs_initial_lat'] = np.min(i_lats) if i_valid else 0
        stats['max_obs_initial_lat'] = np.max(i_lats) if i_valid else 0
        stats['ave_obs_final_lat'] = np.mean(f_lats) if f_valid else 0
        stats['min_obs_final_lat'] = np.min(f_lats) if f_valid else 0
        stats['max_obs_final_lat'] = np.max(f_lats) if f_valid else 0

        stats['intial_lat_by_obs_rs'] = intial_lat_by_obs_rs
        stats['intial_lat_by_obs'] = intial_lat_by_obs
        stats['final_lat_by_obs'] = final_lat_by_obs

        if verbose:
            print('------------------------------')
            print('latencies by observation')
            print("%s: \t\t %f"%('ave_obs_initial_lat_rs',stats['ave_obs_initial_lat_rs']))
            print("%s: \t\t %f"%('std_obs_initial_lat_rs',stats['std_obs_initial_lat_rs']))
            print("%s: %f"%('min_obs_initial_lat_rs',stats['min_obs_initial_lat_rs']))
            print("%s: %f"%('max_obs_initial_lat_rs',stats['max_obs_initial_lat_rs']))
            print("%s: \t\t\t %f"%('ave_obs_initial_lat',stats['ave_obs_initial_lat']))
            print("%s: \t\t\t %f"%('std_obs_initial_lat',stats['std_obs_initial_lat']))
            print("%s: %f"%('min_obs_initial_lat',stats['min_obs_initial_lat']))
            print("%s: %f"%('max_obs_initial_lat',stats['max_obs_initial_lat']))
            print("%s: %f"%('ave_obs_final_lat',stats['ave_obs_final_lat']))
            print("%s: %f"%('min_obs_final_lat',stats['min_obs_final_lat']))
            print("%s: %f"%('max_obs_final_lat',stats['max_obs_final_lat']))

            # for obs in final_lat_by_obs.keys ():
            #     i = intial_lat_by_obs.get(obs,99999.0)
            #     f = final_lat_by_obs.get(obs,99999.0)
            #     print("%s: i %f, f %f"%(obs,i,f))

        return stats

    @staticmethod
    def t_conv_func(t_end,t_start,input_type='datetime',output_units='hours'):
        """ converts input time range to a float difference in desired output units """
        if input_type == "datetime":
            diff = (t_end-t_start).total_seconds ()
        elif input_type == "seconds":
            diff =  t_end-t_start
        else:
            raise NotImplementedError

        if output_units == 'seconds':
            return diff
        elif output_units == 'minutes':
            return diff/60
        elif output_units == 'hours':
            return diff/3600
        else:
            raise NotImplementedError

    @staticmethod
    def calc_av_aoi( d_c_mat, start_calc_window, end_calc_window,input_type="datetime",output_units='hours'):
        """ performs AoI calculation on a formatted input matrix
        
        calculates Age of Information (AoI) based on the data in input matrix.  this is the integration process  performed by summing up the triangular sections representing an AoI curve.

        d_c_mat ( delivery creation matrix) is formatted as such:
            col0        col1
        |    0       |     0      |
        |  t_{d,1}   |  t_{c,1}   |
        |  t_{d,2}   |  t_{c,2}   |
        |  t_{d,3}   |  t_{c,3}   |
        |   ...      |   ...      |    
        |  t_{d,T-1} |  t_{c,T-1} |
        |  t_{d,T}   | dont care  |

        where 'd' stands for delivery, 'c' stands for creation, and T  is the number of time points.

        These points are shown notionally below

        |             /|                         
        |      /|    / |                       /.
        |     / |   /  |  /|                  / .
        |    /  |  /   | / |       (etc)     /  .
        |   /   | /    |/  | /                  .
        |  /    |/     .   |/                   .
        | /     .     ..   .                    .
        |/_____..____._.__..____________________.___
        0      ^^    ^ ^  ^^                    ^
              1a,b   2a,b 3a,b                  Ta
             
        1a: t_{c,1}
        1b: t_{d,1}
        2a: t_{c,2}
        2b: t_{d,2}
        3a: t_{c,3}
        3b: t_{d,3}
        ...
        Ta: t_{d,T}

        When calculating AOI, we sum up a series of triangles defined by the creation and delivery times for data. For each of the time points above, at 'b' new data is being delivered to a destination. At 'a', this data was created. Therefore  the data is already as old as the time between 'a' and 'b' when it arrives at the destination. During the time between 't-1 b' and 't b' the data that was delivered at 't-1 b' is aging, without any updates. At 't b' we receive updated data, that again, already has some age. So we go through all the timepoints summing up the triangle specified by 't-1 b' and 't b', with the tip from 't-1 a' to 't-1 b' subtracted (which actually leaves us with a trapezoid).  once we sum up all of these trapezoids, we divide by the total time to get a time-averaged age of data (or AoI).

        TODO:  update this reference when I  actually publish the equation...
        This calculation is performed per equation ? in ?

        :param d_c_mat: [description]
        :type d_c_mat: [type]
        :param start_calc_window: [description]
        :type start_calc_window: [type]
        :param end_calc_window: [description]
        :type end_calc_window: [type]
        :param input_type: [description], defaults to "datetime"
        :type input_type: str, optional
        :param output_units: [description], defaults to 'hours'
        :type output_units: str, optional
        :returns: [description]
        :rtype: {[type]}
        """

        # Now sum up trapezoidal sections of AoI curve (integration)

        conv_func = GPMetrics.t_conv_func

        aoi_summation = 0
        for t in range(1,len(d_c_mat)):
            trap_addition = (conv_func(d_c_mat[t][0],d_c_mat[t-1][1],input_type,output_units)**2 - conv_func(d_c_mat[t-1][0],d_c_mat[t-1][1],input_type,output_units)**2)/2
            aoi_summation += trap_addition

        av_aoi = aoi_summation / conv_func(end_calc_window,start_calc_window,input_type,output_units)
        return av_aoi

    @staticmethod
    def get_aoi_curve(d_c_mat,base_time,input_type="datetime",x_units='minutes',y_units='hours'):
        """get X and Y for plotting AoI
        
        [description]
        :param d_c_mat: [description]
        :type d_c_mat: [type]
        :param base_time: [description]
        :type base_time: [type]
        :param input_type: [description], defaults to "datetime"
        :type input_type: str, optional
        :param x_units: [description], defaults to 'minutes'
        :type x_units: str, optional
        :param y_units: [description], defaults to 'hours'
        :type y_units: str, optional
        """
        
        conv_func = GPMetrics.t_conv_func

        x = []
        y = []
        for indx, row in  enumerate (d_c_mat):
            if indx==0:
                x.append(conv_func(row[0],base_time,input_type,x_units))
                y.append(conv_func(row[0],row[1],input_type,y_units))
            else:
                last_row = d_c_mat[indx-1]
                x.append(conv_func(row[0],base_time,input_type,x_units))
                y.append(conv_func(row[0],last_row[1],input_type,y_units))
                x.append(conv_func(row[0],base_time,input_type,x_units))
                y.append(conv_func(row[0],row[1],input_type,y_units))

        aoi_curve = {
            'x': x,
            'y': y
        }
        return aoi_curve

    @staticmethod
    def get_av_aoi_routing(d_c_mat_targ,start_calc_window,end_calc_window,dlnk_same_time_slop_s,aoi_units='hours',aoi_plot_t_units='minutes'):
        """ preprocess delivery creation matrix and do AoI calculation, with routing
        
        this code first pre-processes the matrix to get rid of superfluous information that would throw off the AoI calculation. This essentially smooths down the data to the saw-like shape expected for an AoI (versus time) curve.  in the preprocessing we progress through delivery (downlink) times, looking for the earliest creation (observation) time for each delivery time. Here we account for the fact that it can take time to deliver data after its creation
        :param d_c_mat_targ:  the delivery creation matrix for a given target (list of lists, each of two elements - delivery, creation time)
        :type d_c_mat_targ: list(list)
        :param start_calc_window: the start of the window for calculating AoI
        :type start_calc_window: datetime or float
        :param end_calc_window: the end of the window for calculating AoI
        :type end_calc_window: datetime or float
        :param dlnk_same_time_slop_s:  the time delta by which delivery times can differ and still be considered the same time
        :type dlnk_same_time_slop_s: float
        :param output_units:  the time output units used for AoI, defaults to 'hours'
        :type output_units: str, optional
        """

        #  this builds in the assumption that AoI starts at zero time zero
        d_c_mat_filt = [[start_calc_window,start_calc_window]]

        current_time = start_calc_window
        last_creation_time = start_calc_window

        for mat_indx, row in enumerate(d_c_mat_targ):
            if row[0] > current_time:
                current_time = row[0]  # update to most recent delivery (downlink) time

                # the check here with last_creation_time is what ensures that creation (observation) times are always increasing - we want this when looking at data routing for the following reason: it's not helpful to hear later about an earlier event than we previously knew about - because that information does not help us make better TIME SENSITIVE decisions
                if (row[1] > last_creation_time):
                    last_creation_time = row[1]
                    #  add data delivery time (the current time) and the creation time of that data
                    d_c_mat_filt.append([current_time,last_creation_time])

            # if this next row has the same delivery time as previous one
            elif (row[0]-current_time).total_seconds() <= dlnk_same_time_slop_s:

                # but if the creation time of this row is BEFORE what's currently in the d_c mat...
                # (have to check this because we only sorted by delivery time before - not assuming that we also sorted by creation time under each distinct delivery time)
                if (row[1] > last_creation_time):
                    last_creation_time = row[1]
                    d_c_mat_filt[-1][1] = last_creation_time  #replace the last creation time, because it turns out we've seen it more recently

        # add on end time - important for getting proper AoI over whole scenario (first point (start_calc_window,start_calc_window) was already added on matrix)
        d_c_mat_filt.append([end_calc_window,end_calc_window])


        avaoi = GPMetrics.calc_av_aoi( d_c_mat_filt, start_calc_window, end_calc_window,input_type="datetime",output_units=aoi_units)
        aoi_curve = GPMetrics.get_aoi_curve(d_c_mat_filt,start_calc_window,input_type="datetime",x_units=aoi_plot_t_units,y_units=aoi_units)

        return avaoi, aoi_curve

    

    @staticmethod
    def get_av_aoi_no_routing(d_c_mat_targ,start_calc_window,end_calc_window,aoi_units='hours',aoi_plot_t_units='minutes'):
        """ preprocess delivery creation matrix and do AoI calculation, without routing
        
        this code first pre-processes the matrix to get rid of superfluous information that would throw off the AoI calculation. This essentially smooths down the data to the saw-like shape expected for an AoI (versus time) curve.  in the preprocessing we progress through delivery (downlink) times, looking for the earliest creation (observation) time for each delivery time. Here we assume that delivery time equals creation time
        :param d_c_mat_targ:  the delivery creation matrix for a given target (list of lists, each of two elements - delivery, creation time)
        :type d_c_mat_targ: list(list)
        :param start_calc_window: the start of the window for calculating AoI
        :type start_calc_window: datetime or float
        :param end_calc_window: the end of the window for calculating AoI
        :type end_calc_window: datetime or float
        :param dlnk_same_time_slop_s:  the time delta by which delivery times can differ and still be considered the same time
        :type dlnk_same_time_slop_s: float
        :param output_units:  the time output units used for AoI, defaults to 'hours'
        :type output_units: str, optional
        """

        #  this builds in the assumption that AoI starts at zero time zero
        d_c_mat_filt = [[start_calc_window,start_calc_window]]

        last_creation_time = start_calc_window

        for mat_indx, row in enumerate(d_c_mat_targ):
            if row[1] > last_creation_time:
                last_creation_time = row[1]
                # right here we're effectively saying that delivery time is the same as creation time
                # (this will cause a cancellation of the second term in the AoI summation equation)
                d_c_mat_filt.append([last_creation_time,last_creation_time])

        # add on end time - important for getting proper AoI over whole scenario (first point (start_calc_window) was already added on matrix)
        d_c_mat_filt.append([end_calc_window,end_calc_window])

        avaoi = self.calc_av_aoi( d_c_mat_filt, start_calc_window, end_calc_window,input_type="datetime",output_units=aoi_units)
        aoi_curve = self.get_aoi_curve(d_c_mat_filt,start_calc_window,input_type="datetime",x_units=aoi_plot_t_units,y_units=aoi_units)

        return avaoi, aoi_curve

    def preprocess_and_get_aoi(self,rts_by_obs,include_routing,dv_option='scheduled_dv'):
        av_aoi_by_targID = {}
        aoi_curves_by_targID = {}

        # First we need to seperate downlink time and creation time of all obs taken for this target. Put these into a matrix for convenient sorting.
        # for each row of dlnk_obs_times_mat[targ_indx]:
        # column 1 is downlink time
        # column 2 is observation time
        dlnk_obs_times_mat = [[] for targ_indx in range(self.num_targ)]

        # start, center, end...whichever we're using for the latency calculation
        time_option = self.latency_params['dlnk']

        for obs_wind,rts in rts_by_obs.items():

            for targ_ID in obs_wind.target_IDs:

                # skip ignored targets
                if targ_ID in  self.targ_id_ignore_list:
                    continue

                targ_indx = self.all_targ_IDs.index(targ_ID)



                if not include_routing:
                    # add row for this observation. Note: there should be no duplicate observations in obs_winds
                    dlnk_obs_times_mat[targ_indx].append([None,obs_wind.start])

                else:
                    #  want to sort these by earliest time so that we favor earlier downlinks
                    rts.sort (key=lambda rt: getattr(rt.get_dlnk(),time_option))

                    # figure out at which data route we meet the minimum DV downlink requirement
                    cum_dv = 0
                    for dr in rts:
                        if dv_option == 'scheduled_dv':
                            cum_dv += dr.scheduled_dv
                        elif dv_option == 'possible_dv':
                            cum_dv += dr.data_vol
                        else:
                            raise NotImplementedError

                        #  if we have reached our minimum required data volume amount...
                        if cum_dv >= self.min_as_route_dv - self.min_as_obs_dv_slop:

                            dlnk_obs_times_mat[targ_indx].append([getattr(dr.get_dlnk(),time_option), obs_wind.start])
                            #  break because we shouldn't count additional down links from the same observation ( they aren't delivering updated information)
                            break


        for targ_indx in range(self.num_targ):
            dlnk_obs_times_mat_targ = dlnk_obs_times_mat[targ_indx]


            if not include_routing:
                dlnk_obs_times_mat_targ.sort(key=lambda row: row[1])  # sort by creation time

                av_aoi,aoi_curve = self.get_av_aoi_no_routing(dlnk_obs_times_mat_targ, self.met_start_dt, self.met_end_dt,aoi_units=self.aoi_units,aoi_plot_t_units=self.aoi_plot_t_units)

            else:
                dlnk_obs_times_mat_targ.sort(key=lambda row: row[0])  # sort by downlink time

                av_aoi,aoi_curve = self.get_av_aoi_routing(dlnk_obs_times_mat_targ,  self.met_start_dt,self.met_end_dt,self.dlnk_same_time_slop_s,aoi_units=self.aoi_units,aoi_plot_t_units=self.aoi_plot_t_units)
            
            targ_ID = self.all_targ_IDs[targ_indx]
            av_aoi_by_targID[targ_ID] = av_aoi
            aoi_curves_by_targID[targ_ID] = aoi_curve

        return av_aoi_by_targID,aoi_curves_by_targID

    def assess_aoi_by_obs_target(self,rs_routes_by_obs,sched_routes,verbose = True):

        sched_rts_by_obs = self.get_routes_by_obs (sched_routes)

        include_routing=True
        av_aoi_by_targID_rs,aoi_curves_by_targID_rs = self.preprocess_and_get_aoi(rs_routes_by_obs,include_routing,dv_option='possible_dv')
        av_aoi_by_targID_sched,aoi_curves_by_targID_sched = self.preprocess_and_get_aoi(sched_rts_by_obs,include_routing,dv_option='scheduled_dv')

        valid = len(av_aoi_by_targID_rs) > 0

        stats =  {}
        av_aoi_vals_rs = list(av_aoi_by_targID_rs.values())
        av_aoi_vals_sched = list(av_aoi_by_targID_sched.values())
        stats['av_av_aoi_rs'] = np.mean(av_aoi_vals_rs) if valid else None
        stats['av_av_aoi_sched'] = np.mean(av_aoi_vals_sched) if valid else None
        stats['std_av_aoi_rs'] = np.std(av_aoi_vals_rs) if valid else None
        stats['std_av_aoi_sched'] = np.std(av_aoi_vals_sched) if valid else None
        stats['min_av_aoi_rs'] = np.min(av_aoi_vals_rs) if valid else None
        stats['min_av_aoi_sched'] = np.min(av_aoi_vals_sched) if valid else None
        stats['max_av_aoi_rs'] = np.max(av_aoi_vals_rs) if valid else None
        stats['max_av_aoi_sched'] = np.max(av_aoi_vals_sched) if valid else None

        stats['av_aoi_by_targID_rs'] = av_aoi_by_targID_rs
        stats['av_aoi_by_targID_sched'] = av_aoi_by_targID_sched
        stats['aoi_curves_by_targID_rs'] = aoi_curves_by_targID_rs
        stats['aoi_curves_by_targID_sched'] = aoi_curves_by_targID_sched

        if verbose:
            print('------------------------------')
            print('AoI values')
            print("%s: \t\t\t\t %f"%('av_av_aoi_rs',stats['av_av_aoi_rs']))
            print("%s: %f"%('std_av_aoi_rs',stats['std_av_aoi_rs']))
            print("%s: %f"%('min_av_aoi_rs',stats['min_av_aoi_rs']))
            print("%s: %f"%('max_av_aoi_rs',stats['max_av_aoi_rs']))
            print("%s: \t\t\t %f"%('av_av_aoi_sched',stats['av_av_aoi_sched']))
            print("%s: %f"%('std_av_aoi_sched',stats['std_av_aoi_sched']))
            print("%s: %f"%('min_av_aoi_sched',stats['min_av_aoi_sched']))
            print("%s: %f"%('max_av_aoi_sched',stats['max_av_aoi_sched']))

            # for targ_ID in av_aoi_by_targID.keys ():
            #     avaoi = av_aoi_by_targID.get(targ_ID,None)
            #     print("targ_ID %d: av aoi %f"%(targ_ID,avaoi))

        return stats

    @staticmethod
    def  get_aoi_results(update_hists, num_entities,aoi_units,t_units):
        av_aoi_vals = []
        av_aoi_by_ent_indx = {}
        aoi_curves_vals = []
        aoi_curves_by_ent_indx = {}

        for ent_indx in range(num_entities):

            update_hist = update_hists[ent_indx]
            d_c_mat = [[t,lut] for t,lut in zip(update_hist.t,update_hist.last_update_time)]

            start_time = d_c_mat[0][0]
            end_time = d_c_mat[-1][0]
            av_aoi = GPMetrics.calc_av_aoi( d_c_mat, start_time, end_time,input_type="seconds",output_units=aoi_units)
            aoi_curve = GPMetrics.get_aoi_curve(d_c_mat,start_time,input_type="seconds",x_units=t_units,y_units=aoi_units)
            
            av_aoi_vals.append(av_aoi)
            aoi_curves_vals.append(aoi_curve)
            av_aoi_by_ent_indx[ent_indx] = av_aoi
            aoi_curves_by_ent_indx[ent_indx] = aoi_curve

        return av_aoi_vals,av_aoi_by_ent_indx,aoi_curves_vals,aoi_curves_by_ent_indx


    def assess_aoi_sat_cmd(self,sats_cmd_update_hist,verbose = True):
        (av_aoi_vals,
            av_aoi_by_sat_indx,
            aoi_curves_vals,
            aoi_curves_by_sat_indx) = self.get_aoi_results(
                sats_cmd_update_hist,
                self.num_sats,
                self.aoi_units,
                self.aoi_plot_t_units)

        valid = len(av_aoi_vals) > 0

        stats =  {}
        stats['av_av_aoi'] = np.mean(av_aoi_vals) if valid else None
        stats['min_av_aoi'] = np.min(av_aoi_vals) if valid else None
        stats['max_av_aoi'] = np.max(av_aoi_vals) if valid else None
        stats['std_av_aoi'] = np.std(av_aoi_vals) if valid else None

        stats['av_aoi_by_sat_indx'] = av_aoi_by_sat_indx
        stats['aoi_curves_by_sat_indx'] = aoi_curves_by_sat_indx

        if verbose:
            print('------------------------------')
            print('Sat CMD AoI values')
            print("%s: \t\t\t\t %f"%('av_av_aoi',stats['av_av_aoi']))
            print("%s: %f"%('min_av_aoi',stats['min_av_aoi']))
            print("%s: %f"%('max_av_aoi',stats['max_av_aoi']))
            print("%s: %f"%('std_av_aoi',stats['std_av_aoi']))

            # for sat_indx in range(self.num_sats):
            #     avaoi = av_aoi_by_sat_indx.get(sat_indx,None)
            #     print("sat_indx %d: av aoi %f"%(sat_indx,avaoi))

        return stats

    def assess_aoi_sat_tlm(self,sats_tlm_update_hist,verbose = True):
        (av_aoi_vals,
            av_aoi_by_sat_indx,
            aoi_curves_vals,
            aoi_curves_by_sat_indx) =  self.get_aoi_results(
                sats_tlm_update_hist,
                self.num_sats,
                self.aoi_units,
                self.aoi_plot_t_units)

        valid = len(av_aoi_vals) > 0

        stats =  {}
        stats['av_av_aoi'] = np.mean(av_aoi_vals) if valid else None
        stats['min_av_aoi'] = np.min(av_aoi_vals) if valid else None
        stats['max_av_aoi'] = np.max(av_aoi_vals) if valid else None
        stats['std_av_aoi'] = np.std(av_aoi_vals) if valid else None

        stats['av_aoi_by_sat_indx'] = av_aoi_by_sat_indx
        stats['aoi_curves_by_sat_indx'] = aoi_curves_by_sat_indx

        if verbose:
            print('------------------------------')
            print('Sat TLM AoI values')
            print("%s: \t\t\t\t %f"%('av_av_aoi',stats['av_av_aoi']))
            print("%s: %f"%('min_av_aoi',stats['min_av_aoi']))
            print("%s: %f"%('max_av_aoi',stats['max_av_aoi']))
            print("%s: %f"%('std_av_aoi',stats['std_av_aoi']))

            # for sat_indx in range(self.num_sats):
            #     avaoi = av_aoi_by_sat_indx.get(sat_indx,None)
            #     print("sat_indx %d: av aoi %f"%(sat_indx,avaoi))

        return stats



    def assess_resource_margin(self,energy_usage,verbose = True):

        e_margin_by_sat_indx = {}
        ave_e_margin = []
        ave_e_margin_prcnt = []
        min_e_margin = []
        min_e_margin_prcnt = []
        max_e_margin = []
        max_e_margin_prcnt = []
        for sat_indx in range (self.num_sats):
            e_margin = [e - self.sats_emin_Wh[sat_indx] for e in energy_usage['e_sats'][sat_indx]]
            max_margin = self.sats_emax_Wh[sat_indx]-self.sats_emin_Wh[sat_indx]
            e_margin_prcnt = [100*(e - self.sats_emin_Wh[sat_indx])/max_margin for e in energy_usage['e_sats'][sat_indx]]

            e_ave = np.mean(e_margin)
            e_ave_prcnt = np.mean(e_margin_prcnt)
            e_max = np.max(e_margin)
            e_max_prcnt = np.max(e_margin_prcnt)
            e_min = np.min(e_margin)
            e_min_prcnt = np.min(e_margin_prcnt)
            e_margin_by_sat_indx[sat_indx] = {
                "ave": e_ave,
                "max": e_max,
                "min": e_min
            }

            ave_e_margin.append(e_ave)
            ave_e_margin_prcnt.append(e_ave_prcnt)
            min_e_margin.append(e_min)
            min_e_margin_prcnt.append(e_min_prcnt)
            max_e_margin.append(e_max)
            max_e_margin_prcnt.append(e_max_prcnt)

        valid = len(ave_e_margin) > 0

        stats =  {}
        stats['av_ave_e_margin'] = np.mean(ave_e_margin) if valid else None
        stats['av_ave_e_margin_prcnt'] = np.mean(ave_e_margin_prcnt) if valid else None
        stats['min_ave_e_margin'] = np.min(ave_e_margin) if valid else None
        stats['min_ave_e_margin_prcnt'] = np.min(ave_e_margin_prcnt) if valid else None
        stats['max_ave_e_margin'] = np.max(ave_e_margin) if valid else None
        stats['max_ave_e_margin_prcnt'] = np.max(ave_e_margin_prcnt) if valid else None
        stats['std_ave_e_margin'] = np.std(ave_e_margin) if valid else None
        stats['std_ave_e_margin_prcnt'] = np.std(ave_e_margin_prcnt) if valid else None
        stats['min_min_e_margin'] = np.min(min_e_margin) if valid else None
        stats['min_min_e_margin_prcnt'] = np.min(min_e_margin_prcnt) if valid else None


        if verbose:
            print('------------------------------')
            print('Sat energy margin values')
            # print("%s: %f"%('av_ave_e_margin',stats['av_ave_e_margin']))
            print("%s: %f%%"%('av_ave_e_margin_prcnt',stats['av_ave_e_margin_prcnt']))
            # print("%s: %f"%('min_ave_e_margin',stats['min_ave_e_margin']))
            print("%s: %f%%"%('min_ave_e_margin_prcnt',stats['min_ave_e_margin_prcnt']))
            # print("%s: %f"%('max_ave_e_margin',stats['max_ave_e_margin']))
            print("%s: %f%%"%('max_ave_e_margin_prcnt',stats['max_ave_e_margin_prcnt']))
            # print("%s: %f"%('std_ave_e_margin',stats['std_ave_e_margin']))
            print("%s: %f%%"%('std_ave_e_margin_prcnt',stats['std_ave_e_margin_prcnt']))
            # print("%s: %f"%('min_min_e_margin',stats['min_min_e_margin']))
            print("%s: %f%%"%('min_min_e_margin_prcnt',stats['min_min_e_margin_prcnt']))

            # for sat_indx in range(self.num_sats):
            #     print("sat_indx %d: av e margin %f"%(sat_indx,e_margin_by_sat_indx[sat_indx]['ave']))
            #     print("sat_indx %d: min e margin %f"%(sat_indx,e_margin_by_sat_indx[sat_indx]['min']))
            #     print("sat_indx %d: max e margin %f"%(sat_indx,e_margin_by_sat_indx[sat_indx]['max']))

        return stats

    def verify_dr_type(self,dr):
        #  data route can be a multiple types. verify is one of those
        if not type(dr) in self.ALLOWED_DR_TYPES:
            raise RuntimeWarning('Data route object is not of expected type. Expected: %s, saw %s'%(self.ALLOWED_DR_TYPES,dr))

    def calc_overlaps(self,routes_by_obs):

        # Note: no longer  implementing the window_overlap_option of "shared_window". the "mutex_window" option is hardcoded below (i.e.  if two routes share the same window at all they are considered overlapping, we don't consider how much data volume is available for that window)

        #  in this case we increment one for each individual route that a given route is found to have an overlap with
        #  this is hardcoded below currently
        if self.overlap_count_option=='single_overlap':
            pass
        #  in this case  we increment by the number of overlapping windows for each individual route that a given route is found to have an overlap with
        elif self.overlap_count_option=='multiple_overlap':
            raise NotImplementedError
        else:
            raise NotImplementedError


        overlap_cnt_by_route = {}
        overlap_rts_by_route = {}

        obs_by_route = {}

        rts_by_xlnk = {}
        rts_by_dlnk = {}
        #  go through all routes and for each window within the route and the route to a dictionary indexed by the window
        for obs_indx, (obs,rts) in enumerate(routes_by_obs.items()): 
            for dr in rts:
                #  verify the data route type to avoid confusion down the road
                self.verify_dr_type(dr)

                overlap_rts_by_route[dr] = []
                obs_by_route.setdefault(dr, obs)
                for wind in dr.get_winds():
                    if type(wind) == XlnkWindow:
                        rts_by_xlnk.setdefault(wind,[]).append(dr)
                    elif type(wind) == DlnkWindow:
                        rts_by_dlnk.setdefault(wind,[]).append(dr)


        #  look through all of the route lists as indexed by cross-link window and wherever we find a list that has more than one data route,  mark an overlap for all of the data routes in that list
        for wind, rts in rts_by_xlnk.items():
            for dr in rts:
                if len(rts) > 1:
                    #  leaving this is a list addition rather than a set addition to reduce memory requirements
                    overlap_rts_by_route[dr] += rts

        #  do the same for down links, if so desired
        if self.overlap_include_dlnks:
            for wind, rts in rts_by_dlnk.items():
                for dr in rts:
                    if len(rts) > 1:
                        overlap_rts_by_route[dr] += rts


        # sanitize this by removing any duplicates,
        for dr,rts in overlap_rts_by_route.items():

            #  for removing overlap counts for routes from the same observation
            if self.overlap_remove_same_obs:
                obs = obs_by_route[dr]
                overlap_rts_by_route[dr] = list(set(rts) - set(routes_by_obs[obs]))
            else:
                #  take the set of the lists in order to get rid of duplicate elements
                overlap_rts_by_route[dr] = list(set(rts))

            #  final overlap count is the number of routes left
            overlap_cnt_by_route[dr] = len(overlap_rts_by_route[dr])

        return overlap_cnt_by_route,overlap_rts_by_route

    def assess_route_overlap(  self,routes_by_obs,verbose=False):
        # overlap_count_option='single_overlap'
        # overlap_count_option='multiple_overlap'
        # overlap_count_option='shared_window'
        # overlap_count_option='mutex_window'

        t_a = time.time()
        overlap_cnt_by_route,overlap_rts_by_route = self.calc_overlaps(routes_by_obs)
        t_b = time.time()
        time_elapsed = t_b-t_a

        non_overlap_rt_cnt_by_obs = {}
        non_overlap_has_xlnk_rt_cnt_by_obs = {}
        for obs_indx,(obs,rts) in  enumerate (routes_by_obs.items()):
            rts_overlap_counts = {dr:overlap_cnt_by_route[dr] for dr in rts}

            non_overlap_rt_cnt_by_obs[obs] = 0
            non_overlap_has_xlnk_rt_cnt_by_obs[obs] = 0
            for dr_indx, dr in  enumerate(rts):
                #  verify the data route type to avoid confusion down the road
                self.verify_dr_type(dr)

                if rts_overlap_counts[dr] == 0:
                    non_overlap_rt_cnt_by_obs[obs] += 1

                    if dr.has_xlnk():
                        non_overlap_has_xlnk_rt_cnt_by_obs[obs] += 1


        rt_cnt = [len(rts) for rts in routes_by_obs.values()]
        non_overlap_rt_cnt = [cnt for cnt in non_overlap_rt_cnt_by_obs.values()]
        non_overlap_has_xlnk_rt_cnt = [cnt for cnt in non_overlap_has_xlnk_rt_cnt_by_obs.values()]

        stats =  {}
        counts = list(overlap_cnt_by_route.values())
        stats['overlap calc time'] = time_elapsed
        stats['total_num_overlaps'] = sum(counts)
        stats['ave_num_overlaps_by_route'] = np.mean(counts)
        stats['ave_num_rts_by_obs'] = np.mean(rt_cnt)
        stats['mdn_num_rts_by_obs'] = np.median(rt_cnt)
        stats['std_num_rts_by_obs'] = np.std(rt_cnt)
        stats['min_num_rts_by_obs'] = np.min(rt_cnt)
        stats['max_num_rts_by_obs'] = np.max(rt_cnt)
        stats['ave_non_overlap_count_by_obs'] = np.mean(non_overlap_rt_cnt)
        stats['ave_non_overlap_count_has_xlnk_by_obs'] = np.mean(non_overlap_has_xlnk_rt_cnt)
        stats['mdn_num_overlaps_by_route'] = np.median(counts)
        stats['mdn_non_overlap_count_by_obs'] = np.median(non_overlap_rt_cnt)
        stats['mdn_non_overlap_count_has_xlnk_by_obs'] = np.median(non_overlap_has_xlnk_rt_cnt)
        stats['std_num_overlaps_by_route'] = np.std(counts)
        stats['std_non_overlap_count_by_obs'] = np.std(non_overlap_rt_cnt)
        stats['std_non_overlap_count_has_xlnk_by_obs'] = np.std(non_overlap_has_xlnk_rt_cnt)
        stats['min_num_overlaps_by_route'] = np.min(counts)
        stats['min_non_overlap_count_by_obs'] = np.min(non_overlap_rt_cnt)
        stats['min_non_overlap_count_has_xlnk_by_obs'] = np.min(non_overlap_has_xlnk_rt_cnt)
        stats['max_num_overlaps_by_route'] = np.max(counts)
        stats['max_non_overlap_count_by_obs'] = np.max(non_overlap_rt_cnt)
        stats['max_non_overlap_count_has_xlnk_by_obs'] = np.max(non_overlap_has_xlnk_rt_cnt)
        stats['num_obs_no_non_overlaps'] = non_overlap_rt_cnt.count(0)
        stats['num_obs_no_non_overlaps_has_xlnk'] = non_overlap_has_xlnk_rt_cnt.count(0)

        if verbose:
            print('------------------------------')
            print('Num obs: %d'%(len(routes_by_obs.keys ())))
            print('Route overlap statistics')
            print("%s: %d"%('total_num_overlaps',stats['total_num_overlaps']))
            print('------ By route')
            print("%s: %f"%('ave_num_overlaps',stats['ave_num_overlaps_by_route']))
            print("%s: %f"%('mdn_num_overlaps',stats['mdn_num_overlaps_by_route']))
            print("%s: %f"%('std_num_overlaps',stats['std_num_overlaps_by_route']))
            print("%s: %f"%('min_num_overlaps',stats['min_num_overlaps_by_route']))
            print("%s: %f"%('max_num_overlaps',stats['max_num_overlaps_by_route']))
            print('------ By obs, num routes')
            print("%s: %f"%('ave_num_rts_by_obs',stats['ave_num_rts_by_obs']))
            print("%s: %f"%('mdn_num_rts_by_obs',stats['mdn_num_rts_by_obs']))
            print("%s: %f"%('std_num_rts_by_obs',stats['std_num_rts_by_obs']))
            print("%s: %f"%('min_num_rts_by_obs',stats['min_num_rts_by_obs']))
            print("%s: %f"%('max_num_rts_by_obs',stats['max_num_rts_by_obs']))
            print('------ By obs, routes with no overlaps')
            print("%s: %d"%('num_obs_no_non_overlaps',stats['num_obs_no_non_overlaps']))
            print("%s: %f"%('ave_non_overlap_count_by_obs',stats['ave_non_overlap_count_by_obs']))
            print("%s: %f"%('mdn_non_overlap_count_by_obs',stats['mdn_non_overlap_count_by_obs']))
            print("%s: %f"%('std_non_overlap_count_by_obs',stats['std_non_overlap_count_by_obs']))
            print("%s: %f"%('min_non_overlap_count_by_obs',stats['min_non_overlap_count_by_obs']))
            print("%s: %f"%('max_non_overlap_count_by_obs',stats['max_non_overlap_count_by_obs']))
            print('------ By obs, routes with no overlaps, route has xlink')
            print("%s: %d"%('num_obs_no_non_overlaps',stats['num_obs_no_non_overlaps_has_xlnk']))
            print("%s: %f"%('ave_non_overlap_count_by_obs',stats['ave_non_overlap_count_has_xlnk_by_obs']))
            print("%s: %f"%('mdn_non_overlap_count_by_obs',stats['mdn_non_overlap_count_has_xlnk_by_obs']))
            print("%s: %f"%('std_non_overlap_count_by_obs',stats['std_non_overlap_count_has_xlnk_by_obs']))
            print("%s: %f"%('min_non_overlap_count_by_obs',stats['min_non_overlap_count_has_xlnk_by_obs']))
            print("%s: %f"%('max_non_overlap_count_by_obs',stats['max_non_overlap_count_has_xlnk_by_obs']))
            # for obs_indx,(obs,rts) in  enumerate (routes_by_obs.items()):
            #     # print("%2d, dv %-7.2f: overlap count %d"%(dr_indx,dr.data_vol,cnt))
            #     print('obs %d: %s'%(obs_indx,obs))
            #     rts.sort(key=lambda dr: overlap_cnt_by_route[dr])
            #     for dr_indx, dr in  enumerate(rts):
            #         # print("   dr %2d overlap count: %d"%(dr_indx,overlap_cnt_by_route[dr]))
            #         # print("   dr %2d overlap count: %d \t %s"%(dr_indx,overlap_cnt_by_route[dr],dr))
            #         print("   dr %2d overlap count: %d \t %s \t %s"%(dr_indx,overlap_cnt_by_route[dr],[drID for drID in overlap_rts_by_route[dr]],dr))

        # from circinus_tools import debug_tools 
        # debug_tools.debug_breakpt()

        return overlap_cnt_by_route,stats


