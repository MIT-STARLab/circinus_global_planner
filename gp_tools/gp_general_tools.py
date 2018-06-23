from circinus_tools.scheduling.custom_window import   ObsWindow,  DlnkWindow, XlnkWindow

PLAN_WIND_OPTIONS = ['whole','fixed']

def wind_in_planning_window(obj,wind,plan_wind_opt='whole'):

    # anywhere in planning window from planning_start_dt to appropriate plan wind end
    if plan_wind_opt == 'whole':
        start_filt_dt = obj.gp_inst_planning_params['planning_start_dt']
    # anywhere in mutable planning window, from planning_fixed_end to appropriate plan wind end
    elif plan_wind_opt == 'mutable':
        start_filt_dt = obj.gp_inst_planning_params['planning_fixed_end_dt']
    else:
        raise RuntimeWarning('plan window option %s not in implemented options %s'%(plan_wind_opt,PLAN_WIND_OPTIONS))

    if type(wind) == ObsWindow:
        return wind.start >= start_filt_dt and wind.end <= obj.gp_inst_planning_params['planning_end_obs_dt']    
    if type(wind) == XlnkWindow:    
        return wind.start >= start_filt_dt and wind.end <= obj.gp_inst_planning_params['planning_end_xlnk_dt']
    if type(wind) == DlnkWindow:    
        return wind.start >= start_filt_dt and wind.end <= obj.gp_inst_planning_params['planning_end_dlnk_dt']

def dr_in_planning_window(obj,dr,plan_wind_opt='whole'):
    # obs_in_window = wind_in_planning_window(obj,obs)

    for wind in dr.get_winds():
        if not wind_in_planning_window(obj,wind,plan_wind_opt):
            return False

    # # note: dlnk planning end time should always be greater than for obs, xlnk!
    # return obs_in_window and dr_end <= obj.gp_inst_planning_params['planning_end_dlnk_dt']

    return True

def filt_routes(obj,rts):
    return [rt for rt in rts if dr_in_planning_window(obj,rt)]


def filt_routes_by_obs(obj,rts_by_obs):
    return {obs:filt_routes(obj,rts) for obs,rts in rts_by_obs.items() if wind_in_planning_window(obj,obs)}
