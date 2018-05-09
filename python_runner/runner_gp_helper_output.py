from circinus_tools.scheduling.custom_window import   ObsWindow,  DlnkWindow, XlnkWindow,  EclipseWindow
from gp_tools.gp_metrics import GPMetrics
from gp_tools.network_sim.gp_network_sim import GPNetSim

def calc_activity_scheduling_results ( gp_runner_inst,obs_winds,dlnk_winds_flat,rs_routes_by_obs,sched_routes, energy_usage):
    gp_met = GPMetrics(gp_runner_inst.params)

    def in_planning_window(wind):
        if type(wind) == ObsWindow:
            return wind.start >= gp_runner_inst.gp_inst_planning_params['planning_start_dt'] and wind.end <= gp_runner_inst.gp_inst_planning_params['planning_end_obs_xlnk_dt']
        if type(wind) == DlnkWindow:
            return wind.start >= gp_runner_inst.gp_inst_planning_params['planning_start_dt'] and wind.end <= gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt']

    num_collectible_obs_winds = sum(1 for winds in obs_winds for obs in winds if in_planning_window(obs))
    total_collectible_DV_all_obs_winds = sum(obs.data_vol for winds in obs_winds for obs in winds  if in_planning_window(obs))
    total_dlnkable_DV_all_dlnk_winds = sum(dlnk.data_vol for winds in dlnk_winds_flat for dlnk in winds if in_planning_window(dlnk))
    rs_output_routes = [rt for rts in rs_routes_by_obs.values() for rt in rts]
    total_throughput_DV_rs_routes = sum(sum(rt.data_vol for rt in rts) for obs, rts in rs_routes_by_obs.items() if in_planning_window(obs))
    total_collectible_DV_rs_routes = sum(min(obs.data_vol,sum(rt.data_vol for rt in rts)) for obs, rts in rs_routes_by_obs.items() if in_planning_window(obs))

    print('------------------------------')
    print('calc_activity_scheduling_results()')
    print('in scheduling window:')
    print('num_collectible_obs_winds')
    print(num_collectible_obs_winds)
    if len(rs_routes_by_obs.keys()) == 0:
        print('no RS routes found')
    else:
        print('len(rs_output_routes)')
        print(len(rs_output_routes))
    print('len(sched_routes)')
    print(len(sched_routes))
    print('total_collectible_DV_all_obs_winds')
    print(total_collectible_DV_all_obs_winds)
    print('total_dlnkable_DV_all_dlnk_winds')
    print(total_dlnkable_DV_all_dlnk_winds)
    if len(rs_routes_by_obs.keys()) > 0:
        print('total_throughput_DV_rs_routes')
        print(total_throughput_DV_rs_routes)
        print('total_collectible_DV_rs_routes')
        print(total_collectible_DV_rs_routes)
    print('weights')
    print(gp_runner_inst.as_params['obj_weights'])
    # dv_stats = gp_met.assess_dv_all_routes (sched_routes,verbose = True)
    dv_obs_stats = gp_met.assess_dv_by_obs (rs_routes_by_obs,sched_routes,verbose = True)
    lat_stats = gp_met.assess_latency_all_routes (sched_routes,verbose = True)
    lat_obs_stats = gp_met.assess_latency_by_obs (rs_routes_by_obs,sched_routes,verbose = True)
    aoi_targ_stats = gp_met.assess_aoi_by_obs_target(rs_routes_by_obs,sched_routes,verbose = True)

    gp_netsim = GPNetSim ( gp_runner_inst.params, gp_runner_inst.io_proc)
    gp_netsim.sim_tlm_cmd_routing(sched_routes, verbose =  False)
    #  this is indexed by sat index
    sats_cmd_update_hist = gp_netsim.get_all_sats_cmd_update_hist()
    aoi_sat_cmd_stats = gp_met.assess_aoi_sat_cmd(sats_cmd_update_hist,verbose = True)
    #  this is  indexed by ground station index
    sats_tlm_update_hist = gp_netsim.get_all_sats_tlm_update_hist()
    aoi_sat_tlm_stats = gp_met.assess_aoi_sat_tlm(sats_tlm_update_hist,verbose = True)

    resource_margin_stats = gp_met.assess_resource_margin(energy_usage,verbose = True)


    plot_outputs = {}
    plot_outputs['rs_targIDs_found'] = aoi_targ_stats['rs_targIDs_found']
    plot_outputs['as_targIDs_found'] = aoi_targ_stats['as_targIDs_found']
    plot_outputs['obs_aoi_curves_by_targID'] = aoi_targ_stats['aoi_curves_by_targID_sched']
    plot_outputs['obs_aoi_curves_by_targID'] = aoi_targ_stats['aoi_curves_by_targID_sched']
    plot_outputs['initial_lat_by_obs_rs'] = lat_obs_stats['initial_lat_by_obs_rs']
    plot_outputs['initial_lat_by_obs'] = lat_obs_stats['initial_lat_by_obs']
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
            gp_runner_inst.gp_plot.plot_winds(
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

    # do a bunch of stuff to extract the windows from all of the sched_routes as indexed by observation
    # note that this stuff is not thewindows from the scheduled sched_routes, but rather the windows from all the route selected in route selection
    #  start
    sel_obs_winds_flat = [set() for  sat_indx  in range  (gp_runner_inst.sat_params['num_sats'])]
    sel_dlnk_winds_flat = [set() for sat_indx  in range (gp_runner_inst.sat_params['num_sats'])]
    sel_xlnk_winds_flat = [set() for sat_indx  in range (gp_runner_inst.sat_params['num_sats'])]

    for rts_indx, (obs,rts) in enumerate (sel_routes_by_obs.items()):
        obs_winds_rt, dlnk_winds_rt, \
        xlnk_winds_rt, _, _ = gp_runner_inst.io_proc.extract_flat_windows (rts)

        for sat_indx in range  (gp_runner_inst.sat_params['num_sats']):
            [sel_obs_winds_flat[sat_indx] .add( wind)  for wind in obs_winds_rt[sat_indx]]
            [sel_dlnk_winds_flat[sat_indx] .add(wind ) for wind in dlnk_winds_rt[sat_indx]]
            [sel_xlnk_winds_flat[sat_indx] .add(wind ) for wind in xlnk_winds_rt[sat_indx]]

    for sat_indx in range  (gp_runner_inst.sat_params['num_sats']):
        sel_obs_winds_flat[sat_indx] = list(sel_obs_winds_flat[sat_indx])
        sel_dlnk_winds_flat[sat_indx] = list(sel_dlnk_winds_flat[sat_indx])
        sel_xlnk_winds_flat[sat_indx] = list(sel_xlnk_winds_flat[sat_indx])
        sel_obs_winds_flat[sat_indx].sort(key=lambda x: x.start)
        sel_dlnk_winds_flat[sat_indx].sort(key=lambda x: x.start)
        sel_xlnk_winds_flat[sat_indx].sort(key=lambda x: x.start)
    # end

    sched_obs_winds_flat, sched_dlnk_winds_flat, \
    sched_xlnk_winds_flat, link_info_by_wind, route_indcs_by_wind = gp_runner_inst.io_proc.extract_flat_windows (sched_routes,copy_windows= False)

    #
    sats_to_include =  [sat_id for sat_id in gp_runner_inst.sat_params['sat_id_order']]
    # sats_to_include =  [sat_id for sat_id in range(20,30)]
    # sats_to_include = [12,13,14,15,16]

    all_obs_winds,all_dlnk_winds_flat,all_xlnk_winds_flat = all_possible_winds

    # plot all winds
    gp_runner_inst.gp_plot.plot_winds(
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
        base_time = gp_runner_inst.scenario_params['start_utc_dt'],
        plot_title = 'All Possible Activities',
        plot_size_inches = (18,12),
        plot_include_labels = gp_runner_inst.plot_params['plot_AS_include_labels'],
        plot_original_times = True,
        show=  False,
        fig_name='plots/test_all_windows.pdf'
    )


    # plot the selected down links and cross-links this
    gp_runner_inst.gp_plot.plot_winds(
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
        base_time = gp_runner_inst.scenario_params['start_utc_dt'],
        plot_title = 'Scheduled Activities',
        plot_size_inches = (18,12),
        plot_include_labels = gp_runner_inst.plot_params['plot_AS_include_labels'],
        plot_original_times = False,
        show=  False,
        fig_name='plots/test_activity_times.pdf'
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

    # # found_targIDs = metrics_plot_inputs['rs_targIDs_found']
    # found_targIDs = metrics_plot_inputs['as_targIDs_found']
    # # targs_to_include = [targ['id'] for targ in gp_runner_inst.obs_params['targets']]
    # # targs_to_include = [0,3,4,7,8]
    # # targs_to_include = range(15)
    # # targs_to_include = found_targIDs[0:10]
    # targs_to_include = found_targIDs

    # gp_runner_inst.gp_plot.plot_obs_aoi(
    #     targs_to_include,
    #     metrics_plot_inputs['obs_aoi_curves_by_targID'],
    #     gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
    #     gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
    #     base_time = gp_runner_inst.scenario_params['start_utc_dt'],
    #     plot_title = 'Observation Target AoI',
    #     plot_size_inches = (18,12),
    #     show=False,
    #     fig_name='plots/test_obs_aoi_plot.pdf'
    # )

    # # plot obs latency histogram
    # gp_runner_inst.gp_plot.plot_histogram(
    #     data=metrics_plot_inputs['initial_lat_by_obs'].values(),
    #     num_bins = 40,
    #     plot_type = 'histogram',
    #     x_title='Latency (mins)',
    #     y_title='Number of observations',
    #     plot_title = 'Histogram of initial latency by obs - scheduled (min dv %.1f Mb)'%(gp_runner_inst.as_params['min_as_route_dv_Mb']), 
    #     plot_size_inches = (12,6),
    #     show=False,
    #     fig_name='plots/obs_lat_sched_hist.pdf'
    # )

    # # plot obs latency histogram
    # gp_runner_inst.gp_plot.plot_histogram(
    #     data=metrics_plot_inputs['initial_lat_by_obs'].values(),
    #     num_bins = 40,
    #     plot_type = 'cdf',
    #     x_title='Latency (mins)',
    #     y_title='Number of observations',
    #     plot_title = 'Histogram of initial latency by obs - scheduled (min dv %.1f Mb)'%(gp_runner_inst.as_params['min_as_route_dv_Mb']), 
    #     plot_size_inches = (12,6),
    #     show=False,
    #     fig_name='plots/obs_lat_sched_cdf.pdf'
    # )

    # # plot obs latency histogram
    # gp_runner_inst.gp_plot.plot_histogram(
    #     data=metrics_plot_inputs['initial_lat_by_obs_rs'].values(),
    #     num_bins = 40,
    #     plot_type = 'histogram',
    #     x_title='Latency (mins)',
    #     y_title='Number of observations',
    #     plot_title = 'Histogram of initial latency by obs - RS output (min dv %.1f Mb)'%(gp_runner_inst.as_params['min_as_route_dv_Mb']), 
    #     plot_size_inches = (12,6),
    #     show=False,
    #     fig_name='plots/obs_lat_rs_hist.pdf'
    # )

    # # plot obs latency histogram
    # gp_runner_inst.gp_plot.plot_histogram(
    #     data=metrics_plot_inputs['initial_lat_by_obs_rs'].values(),
    #     num_bins = 40,
    #     plot_type = 'cdf',
    #     x_title='Latency (mins)',
    #     y_title='Number of observations',
    #     plot_title = 'Histogram of initial latency by obs - RS output (min dv %.1f Mb)'%(gp_runner_inst.as_params['min_as_route_dv_Mb']), 
    #     plot_size_inches = (12,6),
    #     show=False,
    #     fig_name='plots/obs_lat_rs_cdf.pdf'
    # )

    # # sats_to_include =  [sat_p['sat_id'] for sat_p in gp_runner_inst.sat_orbit_params]
    # # sats_to_include =  range(10)
    # aoi_option = 'cmd'
    # gp_runner_inst.gp_plot.plot_sat_tlm_cmd_aoi(
    #     sats_to_include,
    #     metrics_plot_inputs['cmd_aoi_curves_by_sat_indx'],
    #     aoi_option,
    #     gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
    #     gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
    #     base_time = gp_runner_inst.scenario_params['start_utc_dt'],
    #     plot_title = 'Satellite Command Uplink AoI',
    #     plot_size_inches = (18,12),
    #     show=False,
    #     fig_name='plots/test_cmd_aoi_plot.pdf'
    # )

    # aoi_option = 'tlm'
    # gp_runner_inst.gp_plot.plot_sat_tlm_cmd_aoi(
    #     sats_to_include,
    #     metrics_plot_inputs['tlm_aoi_curves_by_sat_indx'],
    #     aoi_option,
    #     gp_runner_inst.gp_inst_planning_params['planning_start_dt'],
    #     gp_runner_inst.gp_inst_planning_params['planning_end_dlnk_dt'],
    #     base_time = gp_runner_inst.scenario_params['start_utc_dt'],
    #     plot_title = 'Satellite Telemetry Downlink AoI',
    #     plot_size_inches = (18,12),
    #     show=False,
    #     fig_name='plots/test_tlm_aoi_plot.pdf'
    # )