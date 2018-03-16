#  contains code for assessing metrics of global planner output
#
# @author Kit Kennedy
#
#  note that a path is the same as a route.

from datetime import datetime
import numpy as np


class GPMetrics():
    """docstring for GPMetrics"""

    def __init__(self, params):
        self.latency_params = params['latency_calculation']
        self.min_obs_dv = params['metrics_min_obs_dv_Mb']

        # the amount by which the minimum data volume is allowed to be lower than self.min_obs_dv
        self.min_obs_dv_slop = self.min_obs_dv*0.01

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

    def assess_dv_by_obs(self, routes, verbose=False):
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

        rts_by_obs = self.get_routes_by_obs (routes)

        dvs_by_obs =  {}
        for obs, rts in rts_by_obs.items ():
            dvs_by_obs[obs] = sum (rt.scheduled_dv for rt in rts)

        dvs = [dv for dv in dvs_by_obs. values ()]
        
        valid = len(dvs) > 0

        stats['total_dv'] = sum(dvs) if valid else 0
        stats['ave_obs_dv'] = np.mean(dvs) if valid else 0
        stats['min_obs_dv'] = np.min(dvs) if valid else 0
        stats['max_obs_dv'] = np.max(dvs) if valid else 0

        stats['dvs_by_obs'] = dvs_by_obs

        if verbose:
            print('data volume by observation')
            print("%s: %f"%('total_dv',stats['total_dv']))
            print("%s: %f"%('ave_obs_dv',stats['ave_obs_dv']))
            print("%s: %f"%('min_obs_dv',stats['min_obs_dv']))
            print("%s: %f"%('max_obs_dv',stats['max_obs_dv']))

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

        if verbose:
            print('latency for routes')
            print("%s: %f"%('ave_lat_mins',stats['ave_lat_mins']))
            print("%s: %f"%('min_lat_mins',stats['min_lat_mins']))
            print("%s: %f"%('max_lat_mins',stats['max_lat_mins']))


        return stats

    def assess_latency_by_obs(self, routes, verbose=False):
        """ assess downlink latency as grouped by observation
        
        less straightforward than latency by route. First we group by observation,  then we find out how long it took to downlink the first minimum desired amount of data for each observation. based on how long this took we determin the latency of downlink for the observation.
        :param routes: [description]
        :type routes: [type]
        :param verbose: [description], defaults to False
        :type verbose: bool, optional
        :returns: [description]
        :rtype: {[type]}
        """

        stats = {}

        rts_by_obs = self.get_routes_by_obs (routes)

        intial_lat_by_obs =  {}
        final_lat_by_obs =  {}
        for obs, rts in rts_by_obs.items ():
            # start, center, end...whichever we're using for the latency calculation
            time_option = self.latency_params['dlnk']

            #  want to sort these by earliest time so that we favor earlier downlinks
            rts.sort (key=lambda rt: getattr(rt.get_dlnk(),time_option))

            #  figure out the latency for the initial minimum DV downlink
            cum_dv = 0
            for dr in rts:
                cum_dv += dr.scheduled_dv
                
                #  if we have reached our minimum required data volume amount to deem the observation downlinked for the purposes of latency calculation...
                if cum_dv >= self.min_obs_dv - self.min_obs_dv_slop :

                    intial_lat_by_obs[obs] = dr.get_latency(
                        'minutes',
                        obs_option = self.latency_params['obs'], 
                        dlnk_option = self.latency_params['dlnk']
                    )

                    #  break so that we don't continue considering the rest of the data volume
                    break

            # figure out the latency for downlink of all observation data that we chose to downlink
            final_lat_by_obs[obs] = rts[-1].get_latency(
                'minutes',
                obs_option = self.latency_params['obs'], 
                dlnk_option = self.latency_params['dlnk']
            )

        i_lats = [lat for lat in intial_lat_by_obs. values ()]
        f_lats = [lat for lat in final_lat_by_obs. values ()]
        
        i_valid = len(i_lats) > 0
        f_valid = len(f_lats) > 0

        stats['ave_obs_initial_lat'] = np.mean(i_lats) if i_valid else 0
        stats['min_obs_initial_lat'] = np.min(i_lats) if i_valid else 0
        stats['max_obs_initial_lat'] = np.max(i_lats) if i_valid else 0
        stats['ave_obs_final_lat'] = np.mean(f_lats) if f_valid else 0
        stats['min_obs_final_lat'] = np.min(f_lats) if f_valid else 0
        stats['max_obs_final_lat'] = np.max(f_lats) if f_valid else 0

        stats['intial_lat_by_obs'] = intial_lat_by_obs
        stats['final_lat_by_obs'] = final_lat_by_obs

        if verbose:
            print('latencies by observation')
            print("%s: %f"%('ave_obs_initial_lat',stats['ave_obs_initial_lat']))
            print("%s: %f"%('min_obs_initial_lat',stats['min_obs_initial_lat']))
            print("%s: %f"%('max_obs_initial_lat',stats['max_obs_initial_lat']))
            print("%s: %f"%('ave_obs_final_lat',stats['ave_obs_final_lat']))
            print("%s: %f"%('min_obs_final_lat',stats['min_obs_final_lat']))
            print("%s: %f"%('max_obs_final_lat',stats['max_obs_final_lat']))

            for obs in final_lat_by_obs.keys ():
                i = intial_lat_by_obs.get(obs,99999.0)
                f = final_lat_by_obs.get(obs,99999.0)
                print("%s: i %f, f %f"%(obs,i,f))

        return stats

    def calc_aoi_integration( d_c_mat, start_calc_window, end_calc_window,input_type="datetime",output_units='hours'):
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

        def conv_func(t_end,t_start):
            """ converts input time range to a float difference in seconds """
            if input_type == "datetime":
                return (t_end-t_start).total_seconds ()
            elif input_type == "seconds":
                return t_end-t_start

        aoi_summation = 0
        for t in range(1,len(d_c_mat)):
            trap_addition = (conv_func(d_c_mat[t][0],d_c_mat[t-1][1]).total_seconds()**2 - conv_func(d_c_mat[t-1][0],d_c_mat[t-1][1])**2)/2
            aoi_summation += trap_addition


        av_aoi = aoi_summation / conv_func(end_calc_window,start_calc_window)
        if output_units=='hours':
            return av_aoi/3600  # in hours
        elif output_units=='seconds':
            return av_aoi   # in seconds
        else:
            raise NotImplementedError

    # def assess_aoi_by_obs_target(routes,num_targ,targ_ignore_list,start_calc_window, end_calc_window, tstep_td,include_routing=False):

    #     # note: adapted from code in comm_constellation_MDO repo, commit e2d82dbace8e43bb81b1b2b69955f0a143bf7c62

    #     target_av_aoi_no_routing = []
    #     target_av_aoi_routing = []

    #     # First we need to seperate downlink time and creation time of all obs taken for this target. Put these into a matrix for convenient sorting.
    #     # for each row of dlnk_obs_times_mat[targ_indx]:
    #     # column 1 is downlink time
    #     # column 2 is observation time
    #     dlnk_obs_times_mat = [[] for targ_indx in range(num_targ)]

    #     rts_by_obs = self.get_routes_by_obs (routes)

    #     # start, center, end...whichever we're using for the latency calculation
    #     time_option = self.latency_params['dlnk']

    #     for obs_wind,rts in rts_by_obs.items():

    #         for targ_ID in obs_wind.target_IDs:

    #             # skip ignored targets
    #             if targ_ID in targ_ignore_list:
    #                 continue

    #             if not include_routing:
    #                 # add row for this observation. Note: there should be no duplicate observations in obs_winds
    #                 dlnk_obs_times_mat[targ_ID].append([None,obs_wind.start])

    #             else:
    #                 #  want to sort these by earliest time so that we favor earlier downlinks
    #                 rts.sort (key=lambda rt: getattr(rt.get_dlnk(),time_option))

    #                 # figure out at which data route we meet the minimum DV downlink requirement
    #                 cum_dv = 0
    #                 for dr in rts:
    #                     cum_dv += dr.scheduled_dv
                        
    #                     #  if we have reached our minimum required data volume amount...
    #                     if cum_dv >= self.min_obs_dv - self.min_obs_dv_slop:

    #                         dlnk_obs_times_mat[targ_ID].append([getattr(dr.get_dlnk(),time_option), obs_wind.start])


    #     FIX FROM HERE
    #     for targ_indx in range(num_targ):
    #         dlnk_obs_times_mat_targ = dlnk_obs_times_mat[targ_indx]

    #         if not include_routing:
    #             dlnk_obs_times_mat_targ.sort(key=lambda row: row[1])  # sort by creation time

    #             av_aoi_no_routing_temp = calcAvAoI_noRouting(dlnk_obs_times_mat_targ,start_calc_window, end_calc_window)
    #             target_av_aoi_no_routing.append(av_aoi_no_routing_temp)

    #         else:
    #             dlnk_obs_times_mat_targ.sort(key=lambda row: row[0])  # sort by downlink time

    #             av_aoi_routing_temp = calcAvAoI_routing(dlnk_obs_times_mat_targ, start_calc_window,end_calc_window,tstep_td)
    #             target_av_aoi_routing.append(av_aoi_routing_temp)


    #     if not include_routing:
    #         avavaoi = np.mean(target_av_aoi_no_routing)
    #         stdavaoi = np.std(target_av_aoi_no_routing)
    #         # print 'target_av_aoi_no_routing: '+str(target_av_aoi_no_routing)
    #         # print 'mean: '+str(avavaoi)

    #         return avavaoi,stdavaoi
    #     else:
    #         avavaoi = np.mean(target_av_aoi_routing)
    #         stdavaoi = np.std(target_av_aoi_routing)
    #         # print 'target_av_aoi_routing: ' + str(target_av_aoi_routing)
    #         # print 'mean: ' + str(np.mean(target_av_aoi_routing))
    #         return avavaoi,stdavaoi
