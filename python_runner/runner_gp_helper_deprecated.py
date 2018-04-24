def run_nominal_route_selection_v1( gp_runner_inst,obs_winds,dlnk_winds_flat,xlnk_winds):
    gp_ps = gprsv1.GPDataRouteSelection ( gp_runner_inst.params)

    print ('nominal route selection v1')

    obj_weights = {
        "total_dv": 0.45,
        "num_paths_sel": 0.1,
        "latency_sf": 0.45,
    }

    obs_indx =0
    #  list of all routes, indexed by obs_indx ( important for pickling)
    all_routes =[]
    all_routes_obs =[]
    all_stats =[]
    route_times_s =[]

    # this should be unique across all data routes
    dr_uid = 0

    for sat_indx in range( gp_runner_inst.sat_params['num_sats']):
        for  index, obs in  enumerate ( obs_winds[sat_indx]):

            print ("sat_indx")
            print (sat_indx)
            print ("obs")
            print ( index)

            # routes,obs,stats,time_elapsed,dr_uid = gp_runner_inst.run_route_selection(gp_ps,obs,dlnk_winds_flat,xlnk_winds,obj_weights,dr_uid)

            # run the route selection algorithm
            gp_ps.make_model (obs,dlnk_winds_flat,xlnk_winds, obj_weights, verbose = True)
            stats =gp_ps.get_stats (verbose = True)
            t_a = time.time()
            gp_ps.solve ()
            t_b = time.time()
            gp_ps.print_sol ()
            adjust_opt = gp_runner_inst.rs_general_params['adjust_window_timing_post_selection']
            routes,dr_uid = gp_ps. extract_routes ( dr_uid,copy_windows= True,adjust_window_timing=adjust_opt,verbose  = True)

            time_elapsed = t_b-t_a
            # end algorithm

            print ('obs.data_vol, total dlnk dv, ratio')
            print (obs.data_vol,sum(dr.data_vol for dr in routes),sum(dr.data_vol for dr in routes)/obs.data_vol)
            print ('min latency, ave latency, max latency')
            latencies = [(rt.route[-1].end - rt.route[0].end).total_seconds()/60 for rt in  routes]
            if len(latencies) > 0:
                print (np.min( latencies),np.mean( latencies),np.max( latencies))
            else:
                print('no routes found')
            print ('time_elapsed')
            print (time_elapsed)

            all_routes.append ( routes)
            all_routes_obs.append ( obs)
            all_stats.append ( stats)
            route_times_s.append ( time_elapsed)

            obs_indx +=1

            if obs_indx >= 1:
                break

        if obs_indx >= 1:
            break

    return all_routes,all_routes_obs,all_stats,route_times_s,obs_indx

def  setup_test( gp_runner_inst,obs_winds,dlnk_winds_flat,xlnk_winds):
    obs_winds_sel = [[] for sat_indx in range (gp_runner_inst.sat_params['num_sats'])]
    # obs_winds_sel[13].append ( obs_winds[13][0])
    # obs_winds_sel[19].append ( obs_winds[19][0])
    # obs_winds_sel[19].append ( obs_winds[19][1])
    # obs_winds_sel[20].append ( obs_winds[20][6])
    # obs_winds_sel[26].append ( obs_winds[26][0])

    if False:
        gp_runner_inst.gp_plot.plot_winds(
            range (gp_runner_inst.sat_params['num_sats']),
            # [13,19,20,26],
            obs_winds_sel,
            [],
            dlnk_winds_flat,
            [],
            [],
            {},
            gp_runner_inst.scenario_params['start_utc_dt'],
            gp_runner_inst.scenario_params['end_utc_dt'],
            plot_title = 'all obs, dlnks',
            plot_size_inches = (18,12),
            plot_include_labels = True,
            show=  False,
            fig_name='plots/temp1.pdf'
        )

    return obs_winds_sel


def run_test_route_selection( gp_runner_inst,obs_winds,dlnk_winds_flat,xlnk_winds):
    gp_rs = GPDataRouteSelection ( gp_runner_inst.params)

    total_dv_weights = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    num_paths_sel_weights = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    latency_sf_weights = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]
    weights_tups = zip(total_dv_weights,num_paths_sel_weights,latency_sf_weights)

    print ('test route selection')

    # obs_winds_sel =  gp_runner_inst.setup_test(obs_winds,dlnk_winds_flat,xlnk_winds)
    obs_winds_sel = [[] for sat_indx in range (gp_runner_inst.sat_params['num_sats'])]
    # obs_winds_sel[13].append ( obs_winds[13][0])
    # obs_winds_sel[19].append ( obs_winds[19][0])
    # obs_winds_sel[19].append ( obs_winds[19][1])
    obs_winds_sel[20].append ( obs_winds[20][6])
    # obs_winds_sel[26].append ( obs_winds[26][0])

    obs_indx =0
    #  list of all routes, indexed by obs_indx ( important for pickling)
    all_routes =[]
    all_routes_obs =[]
    all_stats =[]
    route_times_s =[]

    #  flatten into a simple list of observations
    obs_winds_sel_flat = flat_list = [item for sublist in obs_winds_sel for item in sublist]


    for  index, obs in  enumerate ( obs_winds_sel_flat):
        for  weights_index, weights in  enumerate(weights_tups):
            print ("obs")
            print ( index)
            print (  'weights_index')
            print (  weights_index)

            obj_weights = {
                "total_dv":   weights[0],
                "num_paths_sel":   weights[1],
                "latency_sf":   weights[2],
            }

            routes,obs,stats,time_elapsed = gp_runner_inst.run_route_selection(gp_rs,obs,dlnk_winds_flat,xlnk_winds,obj_weights)

            print ('obs.data_vol, total dlnk dv, ratio')
            print (obs.data_vol,sum(dr.data_vol for dr in routes),sum(dr.data_vol for dr in routes)/obs.data_vol)
            print ('min latency, ave latency, max latency')
            latencies = [(rt.route[-1].end - rt.route[0].end).total_seconds()/60 for rt in  routes]
            print (np.min( latencies),np.mean( latencies),np.max( latencies))
            print ('time_elapsed')
            print (time_elapsed)

            all_routes.append ( routes)
            all_routes_obs.append ( obs)
            all_stats.append ( stats)
            route_times_s.append ( time_elapsed)

        obs_indx +=1

    return all_routes,all_routes_obs,all_stats,route_times_s,obs_indx,weights_tups