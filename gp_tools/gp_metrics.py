#  contains code for assessing metrics of global planner output
# 
# @author Kit Kennedy
#
#  note that a path is the same as a route. 

import numpy as np

class GPMetrics():
    """docstring for GPMetrics"""
    
    def __init__(self,params):
        pass

    def assess_routes_dv ( self,routes):
        stats = {}

        dvs = []
        for dr in routes:
            dvs.append(dr.scheduled_dv)

        valid = len(dvs) > 0

        stats['total_dv'] = sum(dvs)  if valid else 0
        stats['ave_rt_dv'] = np.mean(dvs) if valid else 0
        stats['min_rt_dv'] = np.min(dvs) if valid else 0
        stats['max_rt_dv'] = np.max(dvs) if valid else 0

        return stats
