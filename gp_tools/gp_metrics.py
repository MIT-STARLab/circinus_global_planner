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

