import pickle
from datetime import datetime, timedelta

def pickle_rtsel_s1_stuff(gp_runner_inst,routes_by_obs,all_stats,route_times_s,obs_indx,ecl_winds,window_uid):

    pickle_stuff =  {}
    pickle_stuff['routes_by_obs'] = routes_by_obs
    pickle_stuff['all_stats'] = all_stats
    pickle_stuff['route_times_s'] = route_times_s
    pickle_stuff['params'] =  gp_runner_inst.params
    pickle_stuff['obs_indx'] = obs_indx
    pickle_stuff['ecl_winds'] = ecl_winds
    pickle_stuff['window_uid'] = window_uid
    pickle_name ='pickles/rs_s1_%s_oi%d_%s' %( gp_runner_inst.other_params['new_pickle_file_name_pre'],obs_indx,datetime.utcnow().isoformat().replace (':','_'))
    with open('%s.pkl' % ( pickle_name),'wb') as f:
        pickle.dump(  pickle_stuff,f)

def unpickle_rtsel_s1_stuff( gp_runner_inst):
    rs_s1_pickle = gp_runner_inst.other_params['rs_s1_pickle_input']
    # rs_s1_pickle = rs_s1_pickle if rs_s1_pickle else gp_runner_inst.pickle_params['route_selection_step1_pickle']

    p = pickle.load (open ( rs_s1_pickle,'rb'))
    #  TODO:  uncommon this  if want to reload parameters from file
    # gp_runner_inst.params = p['params']

    return p['routes_by_obs'],p['all_stats'],p['route_times_s'],p['obs_indx'],p['ecl_winds'],p['window_uid']

def pickle_rtsel_s2_stuff(gp_runner_inst,sel_routes_by_obs,ecl_winds,window_uid,stats_rs2_pre,stats_rs2_post,latest_dr_uid):

    pickle_stuff =  {}
    pickle_stuff['sel_routes_by_obs'] = sel_routes_by_obs
    pickle_stuff['ecl_winds'] = ecl_winds
    pickle_stuff['window_uid'] = window_uid
    pickle_stuff['stats_rs2_pre'] = stats_rs2_pre
    pickle_stuff['stats_rs2_post'] = stats_rs2_post
    pickle_stuff['params'] =  gp_runner_inst.params
    pickle_stuff['latest_dr_uid'] =  latest_dr_uid
    pickle_name ='pickles/rs_s2_%s_%s' %( gp_runner_inst.other_params['new_pickle_file_name_pre'],datetime.utcnow().isoformat().replace (':','_'))
    with open('%s.pkl' % ( pickle_name),'wb') as f:
        pickle.dump(  pickle_stuff,f)

def unpickle_rtsel_s2_stuff( gp_runner_inst):
    rs_s2_pickle = gp_runner_inst.other_params['rs_s2_pickle_input']
    # rs_s2_pickle = rs_s2_pickle if rs_s2_pickle else gp_runner_inst.pickle_params['route_selection_step2_pickle']

    p = pickle.load (open ( rs_s2_pickle,'rb'))
    #  TODO:  uncommon this  if want to reload parameters from file
    # gp_runner_inst.params = p['params']

    return p['sel_routes_by_obs'],p['ecl_winds'],p['window_uid'],p['stats_rs2_pre'],p['stats_rs2_post'],p['latest_dr_uid']


def pickle_actsc_stuff(gp_runner_inst,routes_by_obs,ecl_winds,scheduled_routes,energy_usage,data_usage,window_uid,latest_dr_uid):

    pickle_stuff =  {}
    pickle_stuff['routes_by_obs'] = routes_by_obs
    pickle_stuff['ecl_winds'] = ecl_winds
    pickle_stuff['scheduled_routes'] = scheduled_routes
    pickle_stuff['energy_usage'] = energy_usage
    pickle_stuff['data_usage'] = data_usage
    pickle_stuff['window_uid'] = window_uid
    pickle_stuff['latest_dr_uid'] = latest_dr_uid
    pickle_name ='pickles/as_%s' %(datetime.utcnow().isoformat().replace (':','_'))
    with open('%s.pkl' % ( pickle_name),'wb') as f:
        pickle.dump(  pickle_stuff,f)


def unpickle_actsc_stuff( gp_runner_inst):
    as_pickle = gp_runner_inst.other_params['as_pickle_input']

    p = pickle.load (open ( as_pickle,'rb'))

    return p['routes_by_obs'],p['ecl_winds'],p['scheduled_routes'],p['energy_usage'],p['data_usage'],p['window_uid'],p['latest_dr_uid']