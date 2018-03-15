#  contains code for assessing metrics of global planner output
#
# @author Kit Kennedy
#
#  note that a path is the same as a route.

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

