# Basic network simulation
#
# @author Kit Kennedy
#

from datetime import datetime
import numpy as np
# import ipdb

from .sim_entities import NetSimSat, NetSimGS, UpdateHistory
from gp_tools.custom_activity_window import DlnkWindow, XlnkWindow
from gp_tools.schedule_objects import Dancecard

class GPNetSim():
    """docstring for GPMetrics"""

    def __init__(self, gp_params, io_proc):
        """initializes based on parameters
        
        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """

        sat_params = gp_params['gp_orbit_prop_params']['sat_params']
        gs_params = gp_params['gp_orbit_prop_params']['gs_params']
        as_params = gp_params['gp_general_params']['activity_scheduling_params']
        gp_inst_params = gp_params['gp_instance_params']['activity_scheduling_params']
        gp_general_other_params = gp_params['gp_general_params']['other_params']

        self.sched_start_utc_dt  = gp_inst_params['start_utc_dt']
        self.sched_end_utc_dt  = gp_inst_params['end_utc_dt']
        self.resource_delta_t_s  =as_params['resource_delta_t_s']
        self.num_sats=sat_params['num_sats']
        self.num_gs=gs_params['num_stations']
        self.gs_id_ignore_list=gp_general_other_params['gs_id_ignore_list']
        self.all_gs_IDs = [stat['id'] for stat in gs_params['stations']]

        self.io_proc = io_proc

    def sim_tlm_cmd_routing(self,routes, verbose = False):


        obs_winds, dlnk_winds, xlnk_winds, dummy, dummy = self.io_proc.extract_flat_windows (routes,copy_windows=True)

        # construct a set of dance cards for every satellite, 
        # each of which keeps track of all of the activities of satellite 
        # can possibly execute at any given time slice delta T. 
        # activity_dancecards = [Dancecard(self.sched_start_utc_dt,self.sched_end_utc_dt,self.resource_delta_t_s) for sat_indx in range (self.num_sats)]
        # for sat_indx in range (self.num_sats): 
        #     # activity_dancecards[sat_indx].add_winds_to_dancecard(obs_winds[sat_indx])
        #     activity_dancecards[sat_indx].add_winds_to_dancecard(dlnk_winds[sat_indx])
        #     activity_dancecards[sat_indx].add_winds_to_dancecard(xlnk_winds[sat_indx])

        merged_activity_dancecard = Dancecard(self.sched_start_utc_dt,self.sched_end_utc_dt,self.resource_delta_t_s)
        for sat_indx in range (self.num_sats): 
            # activity_dancecards[sat_indx].add_winds_to_dancecard(obs_winds[sat_indx])
            merged_activity_dancecard.add_winds_to_dancecard(dlnk_winds[sat_indx])
            merged_activity_dancecard.add_winds_to_dancecard(xlnk_winds[sat_indx])

        # get the time points that we will iterate through to step through the simulation
        # timepoints is the indices, whereas timepoints_s is the time values in seconds
        #  NOTE: we assume the same time system for every satellite
        timepoints_s = merged_activity_dancecard.get_timepoint_values(units='seconds', time_option ='relative_to_start')

        net_sim_sats = []
        net_sim_gs = []
        self.net_sim_sats = net_sim_sats
        self.net_sim_gs = net_sim_gs
        for sat_indx in range (self.num_sats): 
            net_sim_sats.append (NetSimSat(sat_indx,self.num_sats, self.num_gs,self.sched_start_utc_dt,self.sched_end_utc_dt))
        for gs_indx in range (self.num_gs): 
            net_sim_gs.append (NetSimGS(gs_indx,self.num_sats, self.num_gs,self.sched_start_utc_dt,self.sched_end_utc_dt))


        for tp_indx,tp in  enumerate ( timepoints_s):
            if verbose:
                if tp_indx % 100 == 0:
                    print ('GP network sim at timepoint: %f'% (tp))


            # we haven't been able to perform any cross-links or downlinks as of the first time point (t=0), so skip it
            if tp_indx == 0:
                continue

            # these are the activities we performed right before the current time point
            activities = merged_activity_dancecard.get_objects_pre_timepoint_indx(tp_indx)

            # need to update the  update times on each of the entities involved in these activities.  this means that any exchange that happens after this takes account of the latest time that each entity updated its information
            # note that none of this will happen for ignored ground stations/satellites, because they should not be present in any of the activity windows
            for act in activities:
                if type (act) == DlnkWindow:
                    gs = net_sim_gs[act.gs_indx]
                    sat = net_sim_sats[act.sat_indx]
                    gs.refresh_update_time (tp)
                    sat.refresh_update_time (tp)

                elif type (act) == XlnkWindow:
                    sat1 = net_sim_sats[act.sat_indx]
                    sat2 = net_sim_sats[act.xsat_indx]
                    sat1.refresh_update_time (tp)
                    sat2.refresh_update_time (tp)

            # now that everybody has updated their own update times, we can update  telemetry and command
            for act in activities:
                if type (act) == DlnkWindow:
                    gs = net_sim_gs[act.gs_indx]
                    sat = net_sim_sats[act.sat_indx]
                    gs.pull_ttc(tp,sat)
                    sat.pull_ttc(tp,gs)

                elif type (act) == XlnkWindow:
                    sat1 = net_sim_sats[act.sat_indx]
                    sat2 = net_sim_sats[act.xsat_indx]
                    sat1.pull_ttc(tp,sat2)
                    sat2.pull_ttc(tp,sat1)

        
    def get_all_sats_cmd_update_hist(self):

        all_sats_update_hist = []

        for sat_indx in range(self.num_sats):
            sim_sat = self.net_sim_sats[sat_indx]
            cmd_update_hist = sim_sat.get_merged_gs_update_hist(self.num_gs,self.all_gs_IDs,self.gs_id_ignore_list)
            all_sats_update_hist.append (cmd_update_hist)

        # ipdb.set_trace ()

        return all_sats_update_hist









    