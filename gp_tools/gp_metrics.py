#  contains code for assessing metrics of global planner output
#
# @author Kit Kennedy
#

import time
from datetime import datetime, timedelta
import numpy as np

from circinus_tools  import time_tools as tt
from circinus_tools.scheduling.custom_window import   ObsWindow,  DlnkWindow, XlnkWindow
from circinus_tools.scheduling.routing_objects import DataRoute, DataMultiRoute

from circinus_tools import debug_tools

class GPMetrics():
    """docstring for GPMetrics"""

    ALLOWED_DR_TYPES = [DataRoute, DataMultiRoute]

    def __init__(self, gp_params):
        """initializes based on parameters
        
        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """
        scenario_params = gp_params['orbit_prop_params']['scenario_params']
        sat_params = gp_params['orbit_prop_params']['sat_params']
        gp_inst_planning_params = gp_params['gp_instance_params']['planning_params']
        obs_params = gp_params['orbit_prop_params']['obs_params']
        gp_general_other_params = gp_params['gp_general_params']['other_params']
        metrics_params = gp_params['gp_general_params']['metrics_params']
        plot_params = gp_params['gp_general_params']['plot_params']
        as_params = gp_params['gp_general_params']['activity_scheduling_params']

        self.latency_params = gp_params['gp_general_params']['other_params']['latency_calculation']
        # self.scenario_start_dt  = scenario_params['start_utc_dt']
        # these are used for AoI calculation
        # todo: update these times once receding horizon working...
        self.met_obs_start_dt  = gp_inst_planning_params['planning_start_dt']
        self.met_obs_end_dt  = gp_inst_planning_params['planning_end_obs_dt']
        self.num_sats=sat_params['num_sats']
        self.num_targ = obs_params['num_targets']
        self.all_targ_IDs = [targ['id'] for targ in obs_params['targets']]
        self.targ_id_ignore_list = gp_general_other_params['targ_id_ignore_list']

        self.min_obs_dv_dlnk_req = as_params['min_obs_dv_dlnk_req_Mb']
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

        # the amount by which the minimum data volume is allowed to be lower than self.min_obs_dv_dlnk_req
        self.min_obs_dv_dlnk_req_slop = self.min_obs_dv_dlnk_req*0.01

        # if two downlink times are within this number of seconds, then they are counted as being at the same time for the purposes of AoI calculation
        self.dlnk_same_time_slop_s = scenario_params['timestep_s'] - 1

    @staticmethod
    def get_routes_by_obs(routes):
        rts_by_obs = {}
        for dr in routes:
            obs = dr.get_obs()
            if not obs in rts_by_obs.keys ():
                rts_by_obs[obs] = []
                
            rts_by_obs[obs].append (dr)

        return rts_by_obs


    def verify_dr_type(self,dr):
        #  data route can be a multiple types. verify is one of those
        if not type(dr) in self.ALLOWED_DR_TYPES:
            raise RuntimeWarning('Data route object is not of expected type. Expected: %s, saw %s'%(self.ALLOWED_DR_TYPES,dr))

    def calc_overlaps(self,routes_by_obs,verbose=False):

        num_obs = len(routes_by_obs.keys())

        overlap_cnt_by_route = {}
        overlap_rts_by_route = {}

        if verbose:
            print("calc overlaps")

        obs_by_route = {}

        rts_by_xlnk = {}
        rts_by_dlnk = {}
        #  go through all routes and for each window within the route and the route to a dictionary indexed by the window
        for obs_indx, (obs,rts) in enumerate(routes_by_obs.items()): 

            if verbose:
                if obs_indx%20 == 0:
                    print("obs %d/%d"%(obs_indx,num_obs-1))

            for dr in rts:
                #  verify the data route type to avoid confusion down the road
                self.verify_dr_type(dr)

                # overlap_rts_by_route[dr] = []
                overlap_cnt_by_route[dr] = 0
                obs_by_route.setdefault(dr, obs)
                for wind in dr.get_winds():
                    if type(wind) == XlnkWindow:
                        rts_by_xlnk.setdefault(wind,[]).append(dr)
                    elif type(wind) == DlnkWindow:
                        rts_by_dlnk.setdefault(wind,[]).append(dr)

        #  look through all of the route lists as indexed by cross-link window and wherever we find a list that has more than one data route,  mark an overlap for all of the data routes in that list
        len_rts_by_xlnk = len(rts_by_xlnk.keys())
        for indx, (wind, rts) in enumerate(rts_by_xlnk.items()):

            if len(rts) > 1:
                for dr1 in rts:
                    overlap_cnt_by_route[dr1] += len(rts)-1

        #  do the same for downlinks, if so desired
        if self.overlap_include_dlnks:
            for wind, rts in rts_by_dlnk.items():

                if len(rts) > 1:
                    for dr1 in rts:
                        overlap_cnt_by_route[dr1] += len(rts)-1

        if verbose:
            print("finished calc overlaps")

        return overlap_cnt_by_route,overlap_rts_by_route

    def assess_route_overlap(  self,routes_by_obs,verbose=False):

        t_a = time.time()
        overlap_cnt_by_route,overlap_rts_by_route = self.calc_overlaps(routes_by_obs,verbose=False)
        t_b = time.time()
        time_elapsed = t_b-t_a

        if verbose:
            print("overlap calc time: %f"%(time_elapsed))

        non_overlap_rt_cnt_by_obs = {}
        non_overlap_has_xlnk_rt_cnt_by_obs = {}
        num_rts_with_xlnk_by_obs = {}
        for obs_indx,(obs,rts) in  enumerate (routes_by_obs.items()):

            rts_overlap_counts = {dr:overlap_cnt_by_route[dr] for dr in rts}

            non_overlap_rt_cnt_by_obs[obs] = 0
            non_overlap_has_xlnk_rt_cnt_by_obs[obs] = 0
            num_rts_with_xlnk_by_obs[obs] = 0
            for dr_indx, dr in  enumerate(rts):
                #  verify the data route type to avoid confusion down the road
                self.verify_dr_type(dr)

                if dr.has_xlnk():
                    num_rts_with_xlnk_by_obs[obs] += 1

                if rts_overlap_counts[dr] == 0:
                    non_overlap_rt_cnt_by_obs[obs] += 1

                    if dr.has_xlnk():
                        non_overlap_has_xlnk_rt_cnt_by_obs[obs] += 1


        rt_cnt = [len(rts) for rts in routes_by_obs.values()]
        rt_cnt_xlnk = [cnt for cnt in num_rts_with_xlnk_by_obs.values()]
        rt_cnt_by_obs = {obs:len(rts) for rts in routes_by_obs.values()}
        non_overlap_rt_cnt = [cnt for cnt in non_overlap_rt_cnt_by_obs.values()]
        non_overlap_has_xlnk_rt_cnt = [cnt for cnt in non_overlap_has_xlnk_rt_cnt_by_obs.values()]

        def lat_getter(dr):
            return dr.get_latency(
                    'minutes',
                    obs_option = self.latency_params['obs'], 
                    dlnk_option = self.latency_params['dlnk']
                )

        dv_dlnkable_by_obs = {obs:sum(rt.data_vol for rt in rts) for obs,rts in routes_by_obs.items()}
        # using 0 here if there are no routes is arbitrary. Affects the average min lat calc
        min_lat_by_obs = {obs:(min(lat_getter(rt) for rt in rts) if rts else None) for obs,rts in routes_by_obs.items()}

        valid = any(rt_cnt)

        stats =  {}
        stats['rt_cnt_by_obs'] = rt_cnt_by_obs
        stats['non_overlap_rt_cnt_by_obs'] = non_overlap_rt_cnt_by_obs
        stats['non_overlap_has_xlnk_rt_cnt_by_obs'] = non_overlap_has_xlnk_rt_cnt_by_obs
        stats['dv_dlnkable_by_obs'] = dv_dlnkable_by_obs
        stats['min_lat_by_obs'] = min_lat_by_obs
        stats['num_rts_with_xlnk_by_obs'] = num_rts_with_xlnk_by_obs


        counts = list(overlap_cnt_by_route.values())
        stats['overlap calc time'] = time_elapsed
        stats['total_num_overlaps'] = sum(counts)
        stats['ave_num_overlaps_by_route'] = np.mean(counts) if valid else None
        stats['ave_num_rts_by_obs'] = np.mean(rt_cnt) if valid else None
        stats['mdn_num_rts_by_obs'] = np.median(rt_cnt) if valid else None
        stats['std_num_rts_by_obs'] = np.std(rt_cnt) if valid else None
        stats['min_num_rts_by_obs'] = np.min(rt_cnt) if valid else None
        stats['max_num_rts_by_obs'] = np.max(rt_cnt) if valid else None
        stats['num_obs_zero_rts'] = rt_cnt.count(0)

        stats['ave_num_xlnk_rts_by_obs'] = np.mean(rt_cnt_xlnk) if valid else None
        stats['std_num_xlnk_rts_by_obs'] = np.std(rt_cnt_xlnk) if valid else None
        stats['min_num_xlnk_rts_by_obs'] = np.min(rt_cnt_xlnk) if valid else None
        stats['max_num_xlnk_rts_by_obs'] = np.max(rt_cnt_xlnk) if valid else None

        stats['ave_dv_dlnkable_by_obs'] = np.mean([min(obs.data_vol,dv_dlnkable) for obs,dv_dlnkable in dv_dlnkable_by_obs.items()]) if valid else None
        stats['total_dv_dlnkable'] = sum([min(obs.data_vol,dv_dlnkable) for obs,dv_dlnkable in dv_dlnkable_by_obs.items()]) if valid else None
        stats['ave_min_lat_by_obs'] = np.mean([val for val in min_lat_by_obs.values() if not val is None]) if valid else None

        stats['ave_non_overlap_count_by_obs'] = np.mean(non_overlap_rt_cnt) if valid else None
        stats['ave_non_overlap_count_has_xlnk_by_obs'] = np.mean(non_overlap_has_xlnk_rt_cnt) if valid else None
        stats['ave_percent_unoverlapped_routes_by_obs'] = np.mean([non_overlap_rt_cnt[i]/rt_cnt[i] for i in range(len(rt_cnt)) if rt_cnt[i] > 0]) if valid else None
        stats['ave_percent_unoverlapped_routes_has_xlnk_by_obs'] = np.mean([non_overlap_has_xlnk_rt_cnt[i]/rt_cnt[i] for i in range(len(rt_cnt)) if rt_cnt[i] > 0]) if valid else None

        stats['mdn_num_overlaps_by_route'] = np.median(counts) if valid else None
        stats['mdn_non_overlap_count_by_obs'] = np.median(non_overlap_rt_cnt) if valid else None
        stats['mdn_non_overlap_count_has_xlnk_by_obs'] = np.median(non_overlap_has_xlnk_rt_cnt) if valid else None
        stats['std_num_overlaps_by_route'] = np.std(counts) if valid else None
        stats['std_non_overlap_count_by_obs'] = np.std(non_overlap_rt_cnt) if valid else None
        stats['std_non_overlap_count_has_xlnk_by_obs'] = np.std(non_overlap_has_xlnk_rt_cnt) if valid else None
        stats['min_num_overlaps_by_route'] = np.min(counts) if valid else None
        stats['min_non_overlap_count_by_obs'] = np.min(non_overlap_rt_cnt) if valid else None
        stats['min_non_overlap_count_has_xlnk_by_obs'] = np.min(non_overlap_has_xlnk_rt_cnt) if valid else None
        stats['max_num_overlaps_by_route'] = np.max(counts) if valid else None
        stats['max_non_overlap_count_by_obs'] = np.max(non_overlap_rt_cnt) if valid else None
        stats['max_non_overlap_count_has_xlnk_by_obs'] = np.max(non_overlap_has_xlnk_rt_cnt) if valid else None
        stats['num_obs_no_non_overlaps'] = non_overlap_rt_cnt.count(0)
        stats['num_obs_no_non_overlaps_has_xlnk'] = non_overlap_has_xlnk_rt_cnt.count(0)

        if verbose:
            print('------------------------------')
            print('Num obs: %d'%(len(routes_by_obs.keys ())))
            print('Num routes total: %d'%(sum(rt_cnt)))
            print('Route overlap statistics')

            if not valid:
                print('no routes found, no valid statistics to display')
                return overlap_cnt_by_route,stats

            # v1
            print("%s: %d"%('total_num_overlaps',stats['total_num_overlaps']))
            # print('------ By route')
            # print("%s: %f"%('ave_num_overlaps',stats['ave_num_overlaps_by_route']))
            # print("%s: %f"%('mdn_num_overlaps',stats['mdn_num_overlaps_by_route']))
            # print("%s: %f"%('std_num_overlaps',stats['std_num_overlaps_by_route']))
            # print("%s: %f"%('min_num_overlaps',stats['min_num_overlaps_by_route']))
            # print("%s: %f"%('max_num_overlaps',stats['max_num_overlaps_by_route']))
            # print('------ By obs, num routes')
            # print("%s: %f"%('ave_num_rts_by_obs',stats['ave_num_rts_by_obs']))
            # print("%s: %f"%('mdn_num_rts_by_obs',stats['mdn_num_rts_by_obs']))
            # print("%s: %f"%('std_num_rts_by_obs',stats['std_num_rts_by_obs']))
            # print("%s: %f"%('min_num_rts_by_obs',stats['min_num_rts_by_obs']))
            # print("%s: %f"%('max_num_rts_by_obs',stats['max_num_rts_by_obs']))
            # print('------ By obs, routes with no overlaps')
            # print("%s: %d"%('num_obs_no_non_overlaps',stats['num_obs_no_non_overlaps']))
            # print("%s: %f"%('ave_non_overlap_count_by_obs',stats['ave_non_overlap_count_by_obs']))
            # print("%s: %f"%('mdn_non_overlap_count_by_obs',stats['mdn_non_overlap_count_by_obs']))
            # print("%s: %f"%('std_non_overlap_count_by_obs',stats['std_non_overlap_count_by_obs']))
            # print("%s: %f"%('min_non_overlap_count_by_obs',stats['min_non_overlap_count_by_obs']))
            # print("%s: %f"%('max_non_overlap_count_by_obs',stats['max_non_overlap_count_by_obs']))
            # print('------ By obs, routes with no overlaps, route has xlink')
            # print("%s: %d"%('num_obs_no_non_overlaps',stats['num_obs_no_non_overlaps_has_xlnk']))
            # print("%s: %f"%('ave_non_overlap_count_by_obs',stats['ave_non_overlap_count_has_xlnk_by_obs']))
            # print("%s: %f"%('mdn_non_overlap_count_by_obs',stats['mdn_non_overlap_count_has_xlnk_by_obs']))
            # print("%s: %f"%('std_non_overlap_count_by_obs',stats['std_non_overlap_count_has_xlnk_by_obs']))
            # print("%s: %f"%('min_non_overlap_count_by_obs',stats['min_non_overlap_count_has_xlnk_by_obs']))
            # print("%s: %f"%('max_non_overlap_count_by_obs',stats['max_non_overlap_count_has_xlnk_by_obs']))

            # v2
            print("%s: \t\t %0.2f"%('ave_num_rts_by_obs',stats['ave_num_rts_by_obs']))
            # print("%s: \t\t %0.2f"%('mdn_num_rts_by_obs',stats['mdn_num_rts_by_obs']))
            print("%s: \t\t %0.2f"%('std_num_rts_by_obs',stats['std_num_rts_by_obs']))
            print("%s: \t\t %0.2f"%('min_num_rts_by_obs',stats['min_num_rts_by_obs']))
            print("%s: \t\t %0.2f"%('max_num_rts_by_obs',stats['max_num_rts_by_obs']))
            print("%s: \t\t %d"%('num_obs_zero_rts',stats['num_obs_zero_rts']))
            print("%s: \t %0.2f"%('ave_num_xlnk_rts_by_obs',stats['ave_num_xlnk_rts_by_obs']))
            print("%s: \t %0.2f"%('std_num_xlnk_rts_by_obs',stats['std_num_xlnk_rts_by_obs']))
            print("%s: \t %0.2f"%('min_num_xlnk_rts_by_obs',stats['min_num_xlnk_rts_by_obs']))
            print("%s: \t\t %0.2f %%"%('ave_pcnt_unovrlp_by_o',100*stats['ave_percent_unoverlapped_routes_by_obs']))
            print("%s: \t %0.2f %%"%('ave_pcnt_unovrlp_xlnk_by_o',100*stats['ave_percent_unoverlapped_routes_has_xlnk_by_obs']))
            print("%s: \t\t %0.2f"%('total_dv_dlnkable',stats['total_dv_dlnkable']))
            print("%s: \t\t %0.2f"%('ave_min_lat_by_obs',stats['ave_min_lat_by_obs']))

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


