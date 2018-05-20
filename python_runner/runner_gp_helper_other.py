def validate_unique_windows( gp_runner_inst,obs_winds,dlnk_winds_flat,xlnk_winds,ecl_winds):
    all_wind_ids = set()

    for indx in range(len(obs_winds)):
        for wind in obs_winds[indx]:
            if wind.window_ID in all_wind_ids:
                raise Exception('Found a duplicate unique window ID where it should not have been possible')
            all_wind_ids.add(wind.window_ID)

    for indx in range(len(dlnk_winds_flat)):
        for wind in dlnk_winds_flat[indx]:
            if wind.window_ID in all_wind_ids:
                raise Exception('Found a duplicate unique window ID where it should not have been possible')
            all_wind_ids.add(wind.window_ID)

    for indx in range(len(xlnk_winds)):
        for indx_2 in range(len(xlnk_winds[indx])):
            for wind in xlnk_winds[indx][indx_2]:
                if wind.window_ID in all_wind_ids:
                    raise Exception('Found a duplicate unique window ID where it should not have been possible')
                all_wind_ids.add(wind.window_ID)

    for indx in range(len(ecl_winds)):
        for wind in ecl_winds[indx]:
            if wind.window_ID in all_wind_ids:
                raise Exception('Found a duplicate unique window ID where it should not have been possible')
            all_wind_ids.add(wind.window_ID)


# def validate_unique_window_objects(gp_runner_inst,sel_routes_by_obs,existing_route_data):
#     wind_obj_ids_by_wind_id = {}
#     # from circinus_tools import debug_tools
#     # debug_tools.debug_breakpt()

#     for rts in sel_routes_by_obs.values():
#         for rt in rts:
#             for wind in rt.get_winds():
#                 if wind.window_ID in wind_obj_ids_by_wind_id.keys():
#                     if not id(wind) == wind_obj_ids_by_wind_id[wind.window_ID]:
#                         raise RuntimeWarning('saw a duplicate window object')
#                 else:
#                     wind_obj_ids_by_wind_id[wind.window_ID] = id(wind)

#     for rt in existing_route_data.get('existing_routes',[]):
#         for wind in rt.get_winds():
#             if wind.window_ID in wind_obj_ids_by_wind_id.keys():
#                 if not id(wind) == wind_obj_ids_by_wind_id[wind.window_ID]:
#                     raise RuntimeWarning('saw a duplicate window object')
#             else:
#                 wind_obj_ids_by_wind_id[wind.window_ID] = id(wind)
