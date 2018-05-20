from circinus_tools.scheduling.custom_window import   ObsWindow,  DlnkWindow, XlnkWindow

def wind_in_planning_window(obj,wind,start_time_override=None):
    planning_start_dt = obj.gp_inst_planning_params['planning_start_dt']
    if start_time_override is not None:
        planning_start_dt == start_time_override

    if type(wind) == ObsWindow:
        return wind.start >= planning_start_dt and wind.end <= obj.gp_inst_planning_params['planning_end_obs_dt']    
    if type(wind) == XlnkWindow:    
        return wind.start >= planning_start_dt and wind.end <= obj.gp_inst_planning_params['planning_end_xlnk_dt']
    if type(wind) == DlnkWindow:    
        return wind.start >= planning_start_dt and wind.end <= obj.gp_inst_planning_params['planning_end_dlnk_dt']

def dr_in_planning_window(obj,dr,start_time_override = None):
    # obs_in_window = wind_in_planning_window(obj,obs)

    for wind in dr.get_winds():
        if not wind_in_planning_window(obj,wind,start_time_override):
            return False

    # # note: dlnk planning end time should always be greater than for obs, xlnk!
    # return obs_in_window and dr_end <= obj.gp_inst_planning_params['planning_end_dlnk_dt']

    return True

def filt_routes(obj,rts):
    return [rt for rt in rts if dr_in_planning_window(obj,rt)]


def filt_routes_by_obs(obj,rts_by_obs):
    return {obs:filt_routes(obj,rts) for obs,rts in rts_by_obs.items() if wind_in_planning_window(obj,obs)}
