from copy import deepcopy

from circinus_tools.scheduling.custom_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow
# from gp_tools.gp_metrics import GPMetrics
from circinus_tools.metrics.metrics_calcs import MetricsCalcs
from circinus_tools  import io_tools
from circinus_tools.plotting import plot_tools as pltl
from gp_tools.network_sim.gp_network_sim import GPNetSim
import gp_tools.gp_general_tools as gp_gen

from circinus_tools import debug_tools

def get_metrics_params(gp_runner_inst_params):
    metrics_params = {}

    scenario_params = gp_runner_inst_params['orbit_prop_params']['scenario_params']
    sat_params = gp_runner_inst_params['orbit_prop_params']['sat_params']
    obs_params = gp_runner_inst_params['orbit_prop_params']['obs_params']
    # sim_metrics_params = gp_runner_inst_params['const_sim_inst_params']['sim_metrics_params']
    # sim_plot_params = gp_runner_inst_params['const_sim_inst_params']['sim_plot_params']
    as_params = gp_runner_inst_params['gp_general_params']['activity_scheduling_params']
    gp_general_other_params = gp_runner_inst_params['gp_general_params']['other_params']
    gp_metrics_params = gp_runner_inst_params['gp_general_params']['metrics_params']
    plot_params = gp_runner_inst_params['gp_general_params']['plot_params']
    gp_inst_planning_params = gp_runner_inst_params['gp_instance_params']['planning_params']

    # these are used for AoI calculation
    metrics_params['met_obs_start_dt']  = gp_inst_planning_params['planning_start_dt']
    metrics_params['met_obs_end_dt']  = gp_inst_planning_params['planning_end_obs_dt']

    # gp_inst_planning_params = gp_params['gp_instance_params']['planning_params']
    # gp_general_other_params = gp_params['gp_general_params']['other_params']
    # metrics_params = gp_params['gp_general_params']['metrics_params']
    # plot_params = gp_params['gp_general_params']['plot_params']

    # self.scenario_start_dt  = scenario_params['start_utc_dt']
    metrics_params['num_sats']=sat_params['num_sats']
    metrics_params['num_targ'] = obs_params['num_targets']
    metrics_params['all_targ_IDs'] = [targ['id'] for targ in obs_params['targets']]
    metrics_params['min_obs_dv_dlnk_req'] = as_params['min_obs_dv_dlnk_req_Mb']

    metrics_params['latency_calculation_params'] = gp_general_other_params['latency_calculation']
    metrics_params['targ_id_ignore_list'] = gp_general_other_params['targ_id_ignore_list']
    metrics_params['aoi_units'] = gp_metrics_params['aoi_units']
    metrics_params['aoi_plot_t_units']=plot_params['time_units']

    metrics_params['sats_emin_Wh'] = []
    metrics_params['sats_emax_Wh'] = []        
    for p_params in sat_params['power_params_by_sat_id'].values():
        sat_edot_by_mode,sat_batt_storage,power_units,charge_eff,discharge_eff = io_tools.parse_power_consumption_params(p_params)

    # metrics_params['sats_emin_Wh'] = [p_params['battery_storage_Wh']['e_min'][p_params['battery_option']] for p_params in metrics_params['power_params']]
    # metrics_params['sats_emax_Wh'] = [p_params['battery_storage_Wh']['e_max'][p_params['battery_option']] for p_params in metrics_params['power_params']]
        metrics_params['sats_emin_Wh'].append(sat_batt_storage['e_min'])
        metrics_params['sats_emax_Wh'].append(sat_batt_storage['e_max'])

    metrics_params['timestep_s'] = scenario_params['timestep_s']

    return metrics_params

def calc_activity_scheduling_results ( gp_runner_inst,obs_winds,dlnk_winds_flat,rs_routes_by_obs,sched_routes, energy_usage):
    mc = MetricsCalcs(get_metrics_params(gp_runner_inst.params))

    # gp_met = GPMetrics(gp_runner_inst.params)

    rs_routes_by_obs_filt = gp_gen.filt_routes_by_obs(gp_runner_inst,rs_routes_by_obs)
    sched_routes_filt = gp_gen.filt_routes(gp_runner_inst,sched_routes)

    # todo: this no longer true because we can have routes that are partly outside of the planning window. need to make results calculation more robust
    # these should actually be equal, because we shouldn't have tried to schedule any route that should be filtered!
    # assert(len(sched_routes_filt) == len(sched_routes))

    num_collectible_obs_winds = sum(1 for winds in obs_winds for obs in winds if gp_gen.wind_in_planning_window(gp_runner_inst,obs))
    total_collectible_DV_all_obs_winds = sum(obs.data_vol for winds in obs_winds for obs in winds  if gp_gen.wind_in_planning_window(gp_runner_inst,obs))
    total_dlnkable_DV_all_dlnk_winds = sum(dlnk.data_vol for winds in dlnk_winds_flat for dlnk in winds if gp_gen.wind_in_planning_window(gp_runner_inst,dlnk))
    rs_output_routes = [rt for rts in rs_routes_by_obs_filt.values() for rt in rts]
    total_throughput_DV_rs_routes = sum(sum(rt.data_vol for rt in rts) for obs, rts in rs_routes_by_obs_filt.items())
    total_collectible_DV_rs_routes = sum(min(obs.data_vol,sum(rt.data_vol for rt in rts)) for obs, rts in rs_routes_by_obs_filt.items())

    print('------------------------------')
    print('calc_activity_scheduling_results()')
    print('in scheduling window:')
    print('num_collectible_obs_winds')
    print(num_collectible_obs_winds)
    if len(rs_routes_by_obs_filt.keys()) == 0:
        print('no RS routes found')
    else:
        print('len(rs_output_routes)')
        print(len(rs_output_routes))
    print('len(sched_routes_filt)')
    print(len(sched_routes_filt))
    print('total_collectible_DV_all_obs_winds')
    print(total_collectible_DV_all_obs_winds)
    print('total_dlnkable_DV_all_dlnk_winds')
    print(total_dlnkable_DV_all_dlnk_winds)
    if len(rs_routes_by_obs_filt.keys()) > 0:
        print('total_throughput_DV_rs_routes')
        print(total_throughput_DV_rs_routes)
        print('total_collectible_DV_rs_routes')
        print(total_collectible_DV_rs_routes)
    print('weights')
    print(gp_runner_inst.as_params['obj_weights'])


    time_units = gp_runner_inst.params['gp_general_params']['plot_params']['time_units']

    # dv_stats = gp_met.assess_dv_all_routes (sched_routes_filt,verbose = True)
    dv_obs_stats = mc.assess_dv_by_obs (rs_output_routes,sched_routes_filt,verbose = True)
    # lat_stats = mc.assess_latency_all_routes (sched_routes_filt,verbose = True)
    lat_obs_stats = mc.assess_latency_by_obs (rs_output_routes,sched_routes_filt,verbose = True)
    aoi_targ_stats = mc.assess_aoi_by_obs_target(rs_output_routes,sched_routes_filt,aoi_x_axis_units=time_units,verbose = True)

    gp_netsim = GPNetSim ( gp_runner_inst.params, gp_runner_inst.io_proc)
    gp_netsim.sim_tlm_cmd_routing(sched_routes_filt, verbose =  False)
    #  this is indexed by sat index
    sats_cmd_update_hist = gp_netsim.get_all_sats_cmd_update_hist()
    aoi_sat_cmd_stats = mc.assess_aoi_sat_ttc_option(sats_cmd_update_hist,ttc_option='cmd',aoi_x_axis_units=time_units,verbose = True)
    #  this is  indexed by ground station index
    sats_tlm_update_hist = gp_netsim.get_all_sats_tlm_update_hist()
    aoi_sat_tlm_stats = mc.assess_aoi_sat_ttc_option(sats_cmd_update_hist,ttc_option='tlm',aoi_x_axis_units=time_units,verbose = True)

    resource_margin_stats = mc.assess_resource_margin(energy_usage,verbose = True)


    plot_outputs = {}
    plot_outputs['rs_targIDs_found'] = aoi_targ_stats['poss_targIDs_found']
    plot_outputs['sched_targIDs_found'] = aoi_targ_stats['exec_targIDs_found']
    plot_outputs['obs_aoi_curves_by_targID'] = aoi_targ_stats['aoi_curves_by_targID_exec']
    plot_outputs['initial_lat_by_obs_rs'] = lat_obs_stats['possible_initial_lat_by_obs_exec']
    plot_outputs['initial_lat_by_obs_sched'] = lat_obs_stats['executed_initial_lat_by_obs_exec']
    plot_outputs['cmd_aoi_curves_by_sat_indx'] = aoi_sat_cmd_stats['aoi_curves_by_sat_indx']
    plot_outputs['tlm_aoi_curves_by_sat_indx'] = aoi_sat_tlm_stats['aoi_curves_by_sat_indx']
    return plot_outputs

def pareto_plot(all_routes):
    total_dv_weights = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    num_paths_sel_weights = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    latency_sf_weights = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]
    weights_tups = zip(total_dv_weights,num_paths_sel_weights,latency_sf_weights)
    gp_runner_inst.gp_plot.plot_route_latdv_pareto(all_routes,weights_tups,'plots/obs_winds_20_6_pareto.pdf')

def plot_route_selection_results ( gp_runner_inst,routes_by_obs,dlnk_winds_flat,xlnk_winds_flat,num_obs_to_plot):

    first_dlnk_indx = 10
    num_dlnks_to_plot = 20
    for obs_indx, (obs,routes) in enumerate (routes_by_obs.items()):
        if obs_indx >= num_obs_to_plot:
            break

        print('Plotting obs index %d'%(obs_indx))

        rts_by_dlnk = OrderedDict()

        for dr in routes:
            dlnk = dr.get_dlnk()
            rts_by_dlnk.setdefault(dlnk,[]).append(dr)

        for dlnk_indx, dlnk in enumerate (rts_by_dlnk.keys()):
            if dlnk_indx >= num_dlnks_to_plot:
                break
            if dlnk_indx < first_dlnk_indx:
                continue

            # if dlnk_indx != 33:
            #     continue

            rts = rts_by_dlnk[dlnk]

            # TODO:  this stuff needs to be  changed
            # (sel_obs_winds_flat, sel_dlnk_winds_flat, sel_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind) = gp_runner_inst.io_proc.extract_flat_windows (routes)
            (sel_obs_winds_flat, sel_dlnk_winds_flat, sel_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind) = gp_runner_inst.io_proc.extract_flat_windows (rts)

            obs = rts[0].get_obs()

            sats_to_include =  range (gp_runner_inst.sat_params['num_sats'])
            # sats_to_include = [1,2,3]

            plot_len = max(gp_runner_inst.rs_general_params['wind_filter_duration_s'],gp_runner_inst.rs_general_params['wind_filter_duration_obs_sat_s'])

            #  plot the selected down links and cross-links
            gp_runner_inst.gp_plot.plot_all_sats_acts(
                sats_to_include,
                sel_obs_winds_flat,
                sel_obs_winds_flat,
                dlnk_winds_flat,
                sel_dlnk_winds_flat,
                # [],
                xlnk_winds_flat,
                sel_xlnk_winds_flat,
                # [],
                route_indcs_by_wind,
                gp_runner_inst.scenario_params['start_utc_dt'],
                # obs.start,
                obs.start + timedelta( seconds= plot_len),
                # gp_runner_inst.scenario_params['start_utc_dt'],
                # gp_runner_inst.scenario_params['start_utc_dt'] + timedelta( seconds= gp_runner_inst.rs_general_params['wind_filter_duration_s']),
                # gp_runner_inst.scenario_params['planning_end_dlnk_dt']-timedelta(minutes=200),
                plot_title = 'Route Plot',
                plot_size_inches = (18,12),
                plot_include_labels = gp_runner_inst.plot_params['plot_RS_include_labels'],
                show= False,
                fig_name='plots/test_rs_o{0}_d{1}_6sat.pdf'.format (obs_indx,dlnk_indx)
            )


def  plot_activity_scheduling_results ( gp_runner_inst,all_possible_winds,sel_routes_by_obs,sched_routes,energy_usage,data_usage,ecl_winds,metrics_plot_inputs):

    # do a bunch of stuff to extract the windows from selected (RS) routes, and update their data volumes for convenient display

    #  make copies of all the objects so we don't screw anything up with our fiddling with the data volumes below
    sel_routes_by_obs_copy = deepcopy(sel_routes_by_obs)
    sel_routes = [dmr for rts in sel_routes_by_obs_copy.values() for dmr in rts]
    sel_obs_winds_flat, sel_dlnk_winds_flat, \
    sel_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind = gp_runner_inst.io_proc.extract_flat_windows (sel_routes,copy_windows= False)
    cum_dv_by_wind = {}
    #  routes and calculate naïve scheduled data volumes for every window
    for obs,rts in sel_routes_by_obs_copy.items():
        obs.scheduled_data_vol = min(sum(dmr.data_vol_for_wind(obs) for dmr in rts),obs.data_vol)

        for dmr in rts:
            for dr in dmr.scheduled_dv_by_dr.keys():
                dmr.scheduled_dv_by_dr[dr] = dr.data_vol

            for wind in dmr.get_winds():
                cum_dv_by_wind.setdefault(wind,0)
                cum_dv_by_wind[wind] += dmr.data_vol_for_wind(wind)
                wind.scheduled_data_vol = min(cum_dv_by_wind[wind],wind.data_vol)
                
    #  go back through and update durations from scheduled data volumes
    for obs,rts in sel_routes_by_obs_copy.items():
        # if obs.scheduled_data_vol == 6500:
        #     from circinus_tools import debug_tools
        #     debug_tools.debug_breakpt()

        for dmr in rts:
            for wind in dmr.get_winds():
                wind.update_duration_from_scheduled_dv()


    # now scheduled obs
    sched_obs_winds_flat, sched_dlnk_winds_flat, \
    sched_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind = gp_runner_inst.io_proc.extract_flat_windows (sched_routes,copy_windows= False)

    #
    sats_to_include =  [sat_id for sat_id in gp_runner_inst.sat_params['sat_id_order']]
    # sats_to_include =  [sat_id for sat_id in range(20,30)]
    # sats_to_include = [12,13,14,15,16]

    all_obs_winds,all_dlnk_winds_flat,all_xlnk_winds_flat = all_possible_winds

    # plot the selected down links and cross-links this
    gp_runner_inst.gp_plot.gp_plot_all_sats_acts(
        sats_to_include,
        sel_obs_winds_flat,
        sched_obs_winds_flat,
        sel_dlnk_winds_flat,
        sched_dlnk_winds_flat,
        sel_xlnk_winds_flat,
        sched_xlnk_winds_flat,
        route_indcs_by_wind,
        gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
        # gp_runner_inst.gp_inst_planning_params['planning_start_dt'] + timedelta( seconds= gp_runner_inst.rs_general_params['wind_filter_duration_s']),
        gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
        base_time_dt = gp_runner_inst.scenario_params['start_utc_dt'],
        plot_title = 'Scheduled Acts / RS Acts',
        plot_size_inches = (18,12),
        plot_include_obs_labels = gp_runner_inst.plot_params['plot_AS_include_obs_labels'],
        plot_include_dlnk_labels = gp_runner_inst.plot_params['plot_AS_include_dlnk_labels'],
        plot_include_xlnk_labels = gp_runner_inst.plot_params['plot_AS_include_xlnk_labels'],
        plot_original_times = False,
        show=  False,
        fig_name='plots/test_sched_windows.pdf'
    )

    # plot all winds
    gp_runner_inst.gp_plot.gp_plot_all_sats_acts(
        sats_to_include,
        all_obs_winds,
        all_obs_winds,
        all_dlnk_winds_flat,
        all_dlnk_winds_flat,
        all_xlnk_winds_flat,
        None,
        None,
        gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
        # gp_runner_inst.gp_inst_planning_params['planning_start_dt'] + timedelta( seconds= gp_runner_inst.rs_general_params['wind_filter_duration_s']),
        gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
        base_time_dt = gp_runner_inst.scenario_params['start_utc_dt'],
        plot_title = 'All Acts',
        plot_size_inches = (18,12),
        plot_include_obs_labels = gp_runner_inst.plot_params['plot_AS_include_obs_labels'],
        plot_include_dlnk_labels = gp_runner_inst.plot_params['plot_AS_include_dlnk_labels'],
        plot_include_xlnk_labels = gp_runner_inst.plot_params['plot_AS_include_xlnk_labels'],
        plot_original_times = True,
        show=  False,
        fig_name='plots/test_all_windows.pdf'
    )

    if len(sel_routes_by_obs.keys()) > 0:
        # plot RS winds
        gp_runner_inst.gp_plot.gp_plot_all_sats_acts(
            sats_to_include,
            all_obs_winds,
            sel_obs_winds_flat,
            all_dlnk_winds_flat,
            sel_dlnk_winds_flat,
            all_xlnk_winds_flat,
            sel_xlnk_winds_flat,
            None,
            gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
            # gp_runner_inst.gp_inst_planning_params['planning_start_dt'] + timedelta( seconds= gp_runner_inst.rs_general_params['wind_filter_duration_s']),
            gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
            base_time_dt = gp_runner_inst.scenario_params['start_utc_dt'],
            plot_title = 'RS Acts / All Acts',
            plot_size_inches = (18,12),
            plot_include_obs_labels = gp_runner_inst.plot_params['plot_RS_include_obs_labels'],
            plot_include_dlnk_labels = gp_runner_inst.plot_params['plot_RS_include_dlnk_labels'],
            plot_include_xlnk_labels = gp_runner_inst.plot_params['plot_RS_include_xlnk_labels'],
            plot_original_times = False,
            show=  False,
            fig_name='plots/test_rs_windows.pdf'
        )



    # # gp_runner_inst.gp_plot.plot_data_circles(
    # #     sats_to_include,
    # #     sched_obs_winds_flat,
    # #     sched_obs_winds_flat,
    # #     sched_dlnk_winds_flat,
    # #     sched_dlnk_winds_flat,
    # #     sched_xlnk_winds_flat,
    # #     sched_xlnk_winds_flat,
    # #     route_indcs_by_wind,
    # #     gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
    # #     # gp_runner_inst.gp_inst_planning_params['planning_start_dt'] + timedelta( seconds= gp_runner_inst.rs_general_params['wind_filter_duration_s']),
    # #     gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
    #       # base_time = gp_runner_inst.scenario_params['start_utc_dt'],
    # #     plot_title = 'Activity Data Volumes',
    # #     plot_size_inches = (18,12),
    # #     plot_include_labels = gp_runner_inst.as_params['plot_include_labels'],
    # #     show=  False,
    # #     fig_name='plots/test_data_volume.pdf'
    # # )

    gp_runner_inst.gp_plot.plot_energy_usage(
        sats_to_include,
        energy_usage,
        ecl_winds,
        gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
        gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
        base_time = gp_runner_inst.scenario_params['start_utc_dt'],
        plot_title = 'Energy Storage Utilization',
        plot_size_inches = (18,12),
        show=  False,
        fig_name='plots/test_energy.pdf'
    )

    # note that little blips upward that appear in data storage plot for "just passing through" crosslink pairs (looks like _____^_____ . Isn't that some great asciiart?) can be caused by a slight overlap in timepoints at beg/end of two back to back activities. This should not cause any problems
    gp_runner_inst.gp_plot.plot_data_usage(
        sats_to_include,
        data_usage,
        ecl_winds,
        gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
        gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
        base_time = gp_runner_inst.scenario_params['start_utc_dt'],
        plot_title = 'Data Storage Utilization',
        plot_size_inches = (18,12),
        show=  False,
        fig_name='plots/test_data.pdf'
    )

    # found_targIDs = metrics_plot_inputs['rs_targIDs_found']
    found_targIDs = metrics_plot_inputs['sched_targIDs_found']
    # targs_to_include = [targ['id'] for targ in gp_runner_inst.obs_params['targets']]
    # targs_to_include = [0,3,4,7,8]
    # targs_to_include = range(15)
    # targs_to_include = found_targIDs[0:10]
    targs_to_include = found_targIDs

    gp_runner_inst.gp_plot.plot_obs_aoi(
        targs_to_include,
        metrics_plot_inputs['obs_aoi_curves_by_targID'],
        gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
        gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
        base_time = gp_runner_inst.scenario_params['start_utc_dt'],
        plot_title = 'Observation Target AoI',
        plot_size_inches = (18,12),
        show=False,
        fig_name='plots/test_obs_aoi_plot.pdf'
    )

    # sats_to_include =  [sat_p['sat_id'] for sat_p in gp_runner_inst.sat_orbit_params]
    # sats_to_include =  range(10)
    curves_by_indx = metrics_plot_inputs['cmd_aoi_curves_by_sat_indx']
    cmd_aoi_curves_by_sat_id = {gp_runner_inst.sat_id_order[sat_indx]:curves for sat_indx,curves in curves_by_indx.items()}
    gp_runner_inst.gp_plot.plot_sat_cmd_aoi(
        sats_to_include,
        cmd_aoi_curves_by_sat_id,
        gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
        gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
        base_time = gp_runner_inst.scenario_params['start_utc_dt'],
        plot_title = 'Satellite Command Uplink AoI',
        plot_size_inches = (18,12),
        show=False,
        fig_name='plots/test_cmd_aoi_plot.pdf'
    )

    curves_by_indx = metrics_plot_inputs['tlm_aoi_curves_by_sat_indx']
    tlm_aoi_curves_by_sat_id = {gp_runner_inst.sat_id_order[sat_indx]:curves for sat_indx,curves in curves_by_indx.items()}
    gp_runner_inst.gp_plot.plot_sat_tlm_aoi(
        sats_to_include,
        tlm_aoi_curves_by_sat_id,
        gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
        gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
        base_time = gp_runner_inst.scenario_params['start_utc_dt'],
        plot_title = 'Satellite Telemetry Downlink AoI',
        plot_size_inches = (18,12),
        show=False,
        fig_name='plots/test_tlm_aoi_plot.pdf'
    )

    # plot obs latency histogram
    pltl.plot_histogram(
        data=metrics_plot_inputs['initial_lat_by_obs_sched'].values(),
        num_bins = 40,
        plot_type = 'histogram',
        x_title='Latency (mins)',
        y_title='Number of observations',
        plot_title = 'Histogram of initial latency by obs - scheduled (min dv %.1f Mb)'%(gp_runner_inst.as_params['min_obs_dv_dlnk_req_Mb']), 
        plot_size_inches = (12,6),
        show=False,
        fig_name='plots/obs_lat_sched_hist.pdf',
        plot_fig_extension = 'pdf' 
    )

    # plot obs latency histogram
    pltl.plot_histogram(
        data=metrics_plot_inputs['initial_lat_by_obs_sched'].values(),
        num_bins = 40,
        plot_type = 'cdf',
        x_title='Latency (mins)',
        y_title='Number of observations',
        plot_title = 'Histogram of initial latency by obs - scheduled (min dv %.1f Mb)'%(gp_runner_inst.as_params['min_obs_dv_dlnk_req_Mb']), 
        plot_size_inches = (12,6),
        show=False,
        fig_name='plots/obs_lat_sched_cdf.pdf',
        plot_fig_extension = 'pdf' 
    )

    # plot obs latency histogram
    pltl.plot_histogram(
        data=metrics_plot_inputs['initial_lat_by_obs_rs'].values(),
        num_bins = 40,
        plot_type = 'histogram',
        x_title='Latency (mins)',
        y_title='Number of observations',
        plot_title = 'Histogram of initial latency by obs - RS output (min dv %.1f Mb)'%(gp_runner_inst.as_params['min_obs_dv_dlnk_req_Mb']), 
        plot_size_inches = (12,6),
        show=False,
        fig_name='plots/obs_lat_rs_hist.pdf',
        plot_fig_extension = 'pdf' 
    )

    # plot obs latency histogram
    pltl.plot_histogram(
        data=metrics_plot_inputs['initial_lat_by_obs_rs'].values(),
        num_bins = 40,
        plot_type = 'cdf',
        x_title='Latency (mins)',
        y_title='Number of observations',
        plot_title = 'Histogram of initial latency by obs - RS output (min dv %.1f Mb)'%(gp_runner_inst.as_params['min_obs_dv_dlnk_req_Mb']), 
        plot_size_inches = (12,6),
        show=False,
        fig_name='plots/obs_lat_rs_cdf.pdf',
        plot_fig_extension = 'pdf' 
    )