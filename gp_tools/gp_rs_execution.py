import time
import gp_tools.gp_route_selection_v2 as gprsv2

from collections import namedtuple

RS_Step1Output = namedtuple('RS_Step1Output','obs_wind routes time_elapsed')

def print_pre_run_msg(obs_wind,obs_indx,sat_indx,sat_obs_index,num_obs):
    print ("")
    print ("----------------------------")
    print ("obs_indx %d/%d"%(obs_indx+1,num_obs))
    print ("sat_indx: %d"%(sat_indx))
    print ("sat_obs_index: %d"%(sat_obs_index))
    print ("obs: %s"%(obs_wind))

def run_rs_step1(gp_rs,obs_wind,dlnk_winds_flat,xlnk_winds,verbose=False):
    t_a = time.time()
    routes = gp_rs.run_step1(obs_wind,dlnk_winds_flat,xlnk_winds,verbose=False)
    t_b = time.time()
    time_elapsed = t_b-t_a

    # should be same as 
    return RS_Step1Output(obs_wind,routes,time_elapsed)

class ParallelRSWorkerWrapper():
    """ wraps a route selection worker in the parallel pool to store its own unique local data for processing"""

    def __init__(self, RS_class, gp_params, dlnk_winds_flat,xlnk_winds, num_obs, verbose = False):
        """ initialize copies of the data for the worker"""
        self.dlnk_winds_flat = dlnk_winds_flat
        self.xlnk_winds = xlnk_winds
        self.verbose = verbose
        self.num_obs = num_obs
        self.gp_rs = RS_class(gp_params)

    def __call__(self, obs_inputs):
        """ do the actual call to route selection"""
        obs_wind = obs_inputs[0]
        obs_indx = obs_inputs[1]
        sat_indx = obs_inputs[2]
        sat_obs_index = obs_inputs[3]

        if self.verbose:
            print_pre_run_msg(obs_wind, obs_indx, sat_indx, sat_obs_index, self.num_obs)

        _,routes,time_elapsed = run_rs_step1(self.gp_rs,obs_wind,self.dlnk_winds_flat,self.xlnk_winds,verbose = False)

        return RS_Step1Output(obs_wind,routes,time_elapsed)
