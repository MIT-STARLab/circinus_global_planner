# Contains functionality for turning input data structures into the 
# objects used by the global planner  scheduling module
# 
# 
# @author Kit Kennedy
#
#  note that a path is the same as a route. 

from  datetime import timedelta

from pyomo import environ  as pe
from pyomo import opt  as po


# for Birdseye use
# from cmdline:
# $ birdseye
# http://localhost:7777/
# in source:
# from birdseye import eye
# @eye
# def yo_mama():
# ...



from .routing_objects import DataRoute

class GPDataRouteSelection():
    """docstring for GP route selection"""
    def __init__(self,params):
        self.num_sats=params['num_sats']
        self.num_paths=params['num_paths']
        self.start_utc_dt  =params['start_utc_dt']
        self.end_utc_dt  =params['end_utc_dt']
        self.M_t_s= 1000000 #  ~11.6 days
        self.M_dv_Mb= 1000000000 #  one petabit
        self.min_path_dv =params['min_path_dv_Mb']
        self.solver_max_runtime =params['solver_max_runtime_s']
        self.wind_filter_duration =  timedelta (seconds =params['wind_filter_duration_s'])

        total_duration =(self.end_utc_dt- self.start_utc_dt).total_seconds ()
        if  total_duration  >self.M_t_s:
            raise Exception ('big M value is too small for %f second scheduling window' % ( total_duration))

    @staticmethod
    def get_crosslink_model_objects(xlink_winds,num_paths,num_sats,start_utc_dt):
        tuple_list =[]
        #  start and end times of cross-links, indexed by xlnk i,j,k tuple
        t_start_dict_by_tuple =  {}
        t_end_dict_by_tuple = {}
        dv_dict_by_tuple = {}

        # loop through all satellite index combinations i,j. make elements for i,j as well as j,i direction
        # explicitly indexing by num_sats just to include a bit of error checking
        for sat_indx in  range (num_sats):
            for xsat_indx in  range (sat_indx+1,num_sats):
                for xlnk_indx in  range (len (xlink_winds[sat_indx][xsat_indx])):
                    # i -> j (i,j,k)
                    tup = (sat_indx, xsat_indx,xlnk_indx)
                    tuple_list.append (tup)
                    # j -> i (j,i,k)
                    tup_sym = (xsat_indx, sat_indx,xlnk_indx)
                    tuple_list.append (tup_sym)

                    wind = xlink_winds[sat_indx][xsat_indx][xlnk_indx]

                    t_start_dict_by_tuple[tup] =  (wind.start - start_utc_dt).total_seconds ()
                    t_end_dict_by_tuple[tup] = (wind.end - start_utc_dt).total_seconds ()
                    dv_dict_by_tuple[tup] = wind.data_vol
                    t_start_dict_by_tuple[tup_sym] =  (wind.start - start_utc_dt).total_seconds ()
                    t_end_dict_by_tuple[tup_sym] = (wind.end - start_utc_dt).total_seconds ()
                    dv_dict_by_tuple[tup_sym] = wind.data_vol


        #  now loop through and include paths as an index
        tuple_list_paths = []

        for path_indx in  range (num_paths):
            for sat_indx in  range (num_sats):
                for xsat_indx in  range (sat_indx+1,num_sats):
                    for xlnk_indx in  range (len (xlink_winds[sat_indx][xsat_indx])): 
                        # i -> j (p,i,j,k)
                        tuple_list_paths.append ((path_indx,sat_indx, xsat_indx,xlnk_indx))
                        # j -> i (p,j,i,k)
                        tuple_list_paths.append ((path_indx,xsat_indx, sat_indx,xlnk_indx))

        return tuple_list, tuple_list_paths, t_start_dict_by_tuple, t_end_dict_by_tuple, dv_dict_by_tuple

    @staticmethod
    def get_downlink_model_objects(dlink_winds_flat,num_paths,num_sats,start_utc_dt):
        tuple_list =[]
        #  start and end times of links, indexed by dlnk i,k tuple
        t_start_dict_by_tuple =  {}
        t_end_dict_by_tuple = {}
        dv_dict_by_tuple = {}

        # loop through all satellites and their downlinks
        # explicitly indexing by num_sats just to include a bit of error checking
        for sat_indx in  range (num_sats):

            for dlnk_indx in  range (len (dlink_winds_flat[sat_indx])):
                tup = (sat_indx,dlnk_indx)
                tuple_list.append (tup)
                wind = dlink_winds_flat[sat_indx][dlnk_indx]
                t_start_dict_by_tuple[tup] =  (wind.start - start_utc_dt).total_seconds ()
                t_end_dict_by_tuple[tup] = (wind.end - start_utc_dt).total_seconds ()
                dv_dict_by_tuple[tup] = wind.data_vol

        #  calculate score factors
        score_factors_dict = {}
        inverse_sum = sum (1/t for t in t_end_dict_by_tuple.values())
        for key in t_end_dict_by_tuple.keys ():
            score_factors_dict[key] = 1/t_end_dict_by_tuple[key] / inverse_sum

        largest_factor = max (score_factors_dict.values ())
        # now normalize by the highest factor found
        for key in t_end_dict_by_tuple.keys ():
            score_factors_dict[key] =score_factors_dict[key]/largest_factor

        #  now loop through and include paths as an index
        tuple_list_paths = []

        for path_indx in  range (num_paths):
            for sat_indx in  range (num_sats):

                    for dlnk_indx in  range (len (dlink_winds_flat[sat_indx])):
                        tuple_list_paths.append ((path_indx,sat_indx,dlnk_indx))

        return tuple_list, tuple_list_paths, t_start_dict_by_tuple, t_end_dict_by_tuple, dv_dict_by_tuple, score_factors_dict

    @staticmethod
    def  filter_windows(obs_wind,dlink_winds_flat,xlink_winds,num_sats,end_utc_dt,wind_filter_duration):
        first =  obs_wind.end
        last =  min ( end_utc_dt, first +  wind_filter_duration)

        dlink_winds_flat_filtered = [[] for sat_indx in  range (num_sats)]
        xlink_winds_flat_filtered = [[[] for xsat_indx in  range ( num_sats)] for sat_indx in  range (num_sats)]

        for sat_indx in  range (num_sats):
            for xsat_indx in  range ( num_sats):
                for wind in xlink_winds[sat_indx][xsat_indx]:
                    if  wind.start > first  and  wind.end  <last:
                        xlink_winds_flat_filtered[sat_indx][xsat_indx]. append ( wind)

            for wind in dlink_winds_flat[sat_indx]:
                if  wind.start > first  and  wind.end  <last:
                    dlink_winds_flat_filtered[sat_indx]. append ( wind)

        return dlink_winds_flat_filtered, xlink_winds_flat_filtered

    def make_model ( self,obs_wind,dlink_winds_flat,xlink_winds, verbose = True):
        model = pe.ConcreteModel()

        self.obs_wind = obs_wind
        self.dlink_winds_flat,self.xlink_winds =  self.filter_windows (obs_wind,dlink_winds_flat,xlink_winds, self.num_sats, self.end_utc_dt, self.wind_filter_duration)

        ##############################
        #  Make indices/ subscripts
        ##############################

        # get lists of indices ( as tuples) and start and end times for downlinks
        (dlnk_subscripts, 
            dlnk_path_subscripts,
            dlnk_t_start_dict, 
            dlnk_t_end_dict,
            dlnk_dv_dict,
            dlnk_sfact_dict)  = self.get_downlink_model_objects(
                                        self.dlink_winds_flat,
                                        self.num_paths,
                                        self.num_sats,
                                        self.start_utc_dt)

        # get lists of indices ( as tuples) and start and end times for  crosslinks
        (xlnk_subscripts, 
            xlnk_path_subscripts,
            xlnk_t_start_dict, 
            xlnk_t_end_dict,
            xlnk_dv_dict)  = self.get_crosslink_model_objects(
                                        self.xlink_winds,
                                        self.num_paths,
                                        self.num_sats,
                                        self.start_utc_dt)

        # subscript for all sats i
        model.sats = pe.Set(initialize= range (self.num_sats))
        # subscript for all paths p
        model.paths = pe.Set(initialize= range (self.num_paths))
        #  subscript for each path p passing from satellite i to j through i,j's xlink k
        model.xlnk_path_subscripts = pe.Set(initialize=xlnk_path_subscripts, dimen=4)
        #  subscript for each path p passing from satellite i  to ground through i's dlink k
        model.dlnk_path_subscripts = pe.Set(initialize=dlnk_path_subscripts, dimen=3)
        #  subscript for dlink k  on satellite i
        model.dlnk_subscripts = pe.Set(initialize=dlnk_subscripts, dimen=2)
        #  subscript for xlink k  on satellite i
        model.xlnk_subscripts = pe.Set(initialize=xlnk_subscripts, dimen=3)
        #  subscripts for arrival and departure  times
        model.arrive_depart = pe.Set(initialize= ['a','d'])

        # model.pprint ()'s'
        print (len(xlnk_path_subscripts))
        print (len(dlnk_path_subscripts))
        print (len(dlnk_subscripts))

        ##############################
        #  Make parameters
        ##############################

        #  start from 0  for our time system
        model.par_min_path_dv = pe.Param (initialize=self.min_path_dv)
        model.par_t_start = pe.Param (initialize=0)
        model.par_t_end = pe.Param (initialize= (self.end_utc_dt - self.start_utc_dt).total_seconds ())

        model.par_t_start_dlnk =pe.Param ( model.dlnk_subscripts,initialize =dlnk_t_start_dict)
        model.par_t_end_dlnk =pe.Param ( model.dlnk_subscripts,initialize =dlnk_t_end_dict)
        model.par_dlnk_sf =pe.Param ( model.dlnk_subscripts,initialize =dlnk_sfact_dict)

        model.par_t_start_xlnk =pe.Param ( model.xlnk_subscripts,initialize =xlnk_t_start_dict)
        model.par_t_end_xlnk =pe.Param ( model.xlnk_subscripts,initialize =xlnk_t_end_dict)

        model.par_t_end_obs = pe.Param (initialize=(obs_wind.end - self.start_utc_dt).total_seconds ())
        model.par_obs_occ = pe.Param(model.sats,initialize ={ i: int (obs_wind.sat_indx == i) for i in  model.sats})

        model.par_dlnk_dv = pe.Param (model.dlnk_subscripts,initialize =dlnk_dv_dict)
        model.par_xlnk_dv = pe.Param (model.xlnk_subscripts,initialize =xlnk_dv_dict)

        model.par_obs_dv = pe.Param (initialize = obs_wind.data_vol)
        
        # relative weightings for the objective terms
        model.par_obj_weight1 = pe.Param (initialize = 0.495)
        model.par_obj_weight2 = pe.Param (initialize = 0.01)
        model.par_obj_weight3 = pe.Param (initialize = 0.495)

        #  quick sanity check
        for dv in dlnk_dv_dict.values ():
            if dv >self.M_dv_Mb/1000:
                raise  Exception (' value of self.M_dv_Mb is too small')

        ##############################
        #  Make variables
        ##############################

        #  path occupancy variables  for each crosslink x_(p,i,j,k) [1]
        model.var_xlnk_path_occ  = pe.Var (model.xlnk_path_subscripts, within =pe.Binary)
        #  path occupancy variables  for each downlink d_(p,i,k) [2]
        model.var_dlnk_path_occ  = pe.Var (model.dlnk_path_subscripts, within =pe.Binary)
        #  data volume throughput for each path v_p [3]
        model.var_path_dv  = pe.Var ( model.paths, within = pe.NonNegativeReals)
        #  data volume throughput for path p assigned to downlink d_(i,k) [6]
        model.var_path_dv_dlnk  = pe.Var ( model.paths,  model.dlnk_subscripts,  within = pe.NonNegativeReals)
        #  indicator that path p has been selected I_p [4]
        model.var_path_indic  = pe.Var (model.paths, within =pe.Binary)
        #  arrival and departure times for each path p through sat i  t_(p,i,a), t_(p,i,d) [5]
        model.var_time_a_d   = pe.Var (model.paths,  model.sats,  model.arrive_depart, within =pe.NonNegativeReals)
        #  indicates that no paths could be constructed at all [8]
        model.var_route_fail  = pe.Var( model.sats, within=pe.Binary)
 
        ##############################
        #  Make constraints
        ##############################

        # TODO: renumber  these with the final numbering

        # path p arrival time at sat i  must be greater than t start
        def c1_rule( model,p,i):
            return model.var_time_a_d [p,i,'a']  >= model.par_t_start
        model.c1 =pe.Constraint ( model.paths, model.sats,  rule=c1_rule)

        def c2_rule( model,p,i):
            return model.var_time_a_d [p,i,'d'] <= model.par_t_end
        model.c2 =pe.Constraint ( model.paths, model.sats,  rule=c2_rule)

        def c3_rule( model,p,i):
            return model.var_time_a_d [p,i,'a']  <= model.var_time_a_d [p,i,'d']
        model.c3 =pe.Constraint ( model.paths, model.sats,  rule=c3_rule)

        def c4_rule( model,p,i,j,k):
            return  model.var_time_a_d[p,j,'a']  >= model.var_xlnk_path_occ [p,i,j,k]*model.par_t_end_xlnk[i,j,k]
        model.c4 =pe.Constraint ( model.xlnk_path_subscripts,  rule=c4_rule)

        def c5_rule( model,p,i,j,k):
            return  model.var_time_a_d[p,i,'d']  <= model.par_t_start_xlnk[i,j,k] +  self.M_t_s *(1-model.var_xlnk_path_occ [p,i,j,k])
        model.c5 =pe.Constraint ( model.xlnk_path_subscripts,  rule=c5_rule)

        def c6_rule( model,p,i):
            return  model.var_time_a_d[p,i,'a']  >= model.par_t_end_obs
        model.c6 =pe.Constraint ( model.paths, model.sats,  rule=c6_rule)

        def c7_rule( model,p,i,k):
            return  model.var_time_a_d[p,i,'d']  <= model.par_t_start_dlnk[i,k]  +self.M_t_s *(1-model.var_dlnk_path_occ [p,i,k])
        model.c7 =pe.Constraint ( model.dlnk_path_subscripts,  rule=c7_rule)
        
        def c8b_rule( model,i,j,k):
            return  sum(model.var_xlnk_path_occ[p,i,j,k] for p in  model.paths) <= 1
        model.c8b =pe.Constraint ( model.xlnk_subscripts,  rule=c8b_rule)

        def c11_rule( model,p):
            return sum([model.var_dlnk_path_occ[tuple([p])+dlnk_subsc] for dlnk_subsc in  model.dlnk_subscripts]) <= 1
        model.c11 =pe.Constraint (  model.paths,  rule =c11_rule )

        #  constraint 12
        def c12_rule( model,p):
            return  model.var_path_dv[p]  <= model.par_obs_dv
        model.c12 =pe.Constraint (   model.paths,  rule = c12_rule )

        def c13_rule( model,p,i,j,k):
            return model.var_path_dv[p] <=  model.par_xlnk_dv[i,j,k]  +  self.M_dv_Mb*(1-model.var_xlnk_path_occ[p,i,j,k])
        model.c13 =pe.Constraint (   model.xlnk_path_subscripts,  rule = c13_rule )

        def c15b_rule( model,p):
            return sum (model.var_path_dv_dlnk[tuple([p])+dlnk_subsc] for dlnk_subsc in  model.dlnk_subscripts) >= model.par_min_path_dv*model.var_path_indic[p]
        model.c15b =pe.Constraint ( model.paths, rule=c15b_rule)

        def c16_rule( model,i,k): 
            return sum (model.var_path_dv_dlnk[p,i,k] for p in model.paths) <= model.par_dlnk_dv[i,k]
        model.c16 =pe.Constraint ( model.dlnk_subscripts,rule=c16_rule)

        def c17_rule( model,p): 
            return sum (model.var_path_dv_dlnk[tuple([p])+dlnk_subsc] for dlnk_subsc in  model.dlnk_subscripts)  <=  model.var_path_dv[p]
        model.c17 =pe.Constraint ( model.paths,rule=c17_rule)

        def c18_rule( model,p,i,k): 
            return model.var_path_dv_dlnk[p,i,k]  <=  model.var_dlnk_path_occ[p,i,k]*model.par_dlnk_dv[i,k]
        model.c18 =pe.Constraint ( model.paths, model.dlnk_subscripts,rule=c18_rule)

 
        #  constraint  20
        model.c20  = pe.ConstraintList()
        for p in model.paths: 
            for j in model.sats:
                j_terms_arrive = []
                j_terms_depart = []

                #  could have been a list comprehension but whatever
                for k_d in range(len (self.dlink_winds_flat[j])):
                    j_terms_depart.append (model.var_dlnk_path_occ[p,j,k_d] )

                for i in model.sats:
                    if i == j:
                        continue

                    #  use this to account for the fact that the cross-link Windows matrix is upper triangular
                    num_i_j_xlnks = max(len(self.xlink_winds[i][j]),len(self.xlink_winds[j][i]))

                    j_terms_arrive += [
                        model.var_xlnk_path_occ[p,i,j,k] 
                        for k in range (num_i_j_xlnks)
                    ]
                    j_terms_depart += [
                        model.var_xlnk_path_occ[p,j,i,k] 
                        for k in range (num_i_j_xlnks)
                    ]

                j_terms_arrive.append (model.par_obs_occ[j])

                # add the term accounting for possible failure to route any paths
                j_terms_depart.append (model.var_route_fail[j])

                model.c20.add( sum(j_terms_arrive) == sum(j_terms_depart))


        def c21_rule( model,p):
            return sum([model.var_dlnk_path_occ[tuple([p])+dlnk_subsc] for dlnk_subsc in  model.dlnk_subscripts]) <= model.var_path_indic[p]
        model.c21 =pe.Constraint (  model.paths,  rule =c21_rule )

        def c22_rule( model,p,i):
            return  model.var_route_fail[i] <= 1 - model.var_path_indic[p]
        model.c22 =pe.Constraint (  model.paths,  model.sats, rule =c22_rule )

        # constraint 23
        model.c23  = pe.ConstraintList()
        for j in model.sats:
            if j != self.obs_wind.sat_indx:
                model.c23.add( model.var_route_fail[i] <= 0)

        ##############################
        #  Make objective
        ##############################

        def obj_rule(model):
            return (
                model.par_obj_weight1 * 1/self.num_paths/model.par_obs_dv * sum(model.var_path_dv_dlnk[p,i,k] for i,k in model.dlnk_subscripts for p in model.paths )  +
                model.par_obj_weight2 * 1/self.num_paths                  * sum(model.var_path_indic[p] for p in model.paths) +
                model.par_obj_weight3 * 1/self.num_paths                  * sum(model.var_dlnk_path_occ[p,i,k]*model.par_dlnk_sf[i,k] for i,k in model.dlnk_subscripts for p in model.paths )
            )
        model.obj = pe.Objective( rule=obj_rule, sense=pe.maximize )

        self.model = model

    def get_stats(self,verbose=True):
        stats = {}
        stats['num_dlnk_windows'] = sum([len (self.dlink_winds_flat[sat_indx]) for sat_indx in range (self.num_sats)])
        stats['num_xlnk_windows'] = sum([len (self.xlink_winds[sat_indx][xsat_indx]) for xsat_indx in range (self.num_sats) for sat_indx in range (self.num_sats)])
        stats['num_variables'] = self.model.nvariables ()
        stats['num_nobjectives'] = self.model.nobjectives ()
        stats['num_nconstraints'] = self.model.nconstraints ()

        if verbose:
            print ( "Obs dv: %f" % ( self.obs_wind.data_vol))
            print ( "Considering %d downlink windows" % (stats['num_dlnk_windows']))
            print ( "Considering %d crosslink windows" % (stats['num_xlnk_windows']))
            print ( 'self.model.nvariables ()')
            print ( self.model.nvariables ())
            print ( 'self.model.nobjectives ()')
            print ( self.model.nobjectives ())
            print ( 'self.model.nconstraints ()')
            print ( self.model.nconstraints ())

        return stats




    # lifted from  Jeff Menezes' code at https://github.mit.edu/jmenezes/Satellite-MILP/blob/master/sat_milp_pyomo.py
    def solve(self):

        solver = po.SolverFactory('gurobi')
        results =  solver.solve(self.model, tee=True, keepfiles=False, options_string=" time_limit=%f" % ( self.solver_max_runtime))
        
        # opt = po.SolverFactory("gurobi")
        # solver_manager = po.SolverManagerFactory('neos')
        # results = solver_manager.solve(self.model, opt=opt, options_string=" time_limit=%f" % ( self.solver_max_runtime))

        if (results.solver.status == po.SolverStatus.ok) and (results.solver.termination_condition == po.TerminationCondition.optimal):
            print('this is feasible and optimal')
        elif results.solver.termination_condition == po.TerminationCondition.infeasible:
            print ('infeasible')
        else:
            # something else is wrong
            print (results.solver)

    # lifted from  Jeff Menezes' code at https://github.mit.edu/jmenezes/Satellite-MILP/blob/master/sat_milp_pyomo.py
    def print_sol_all(self):
        for v in self.model.component_objects(pe.Var, active=True):
            print ("Variable",v)
            varobject = getattr(self.model, str(v))
            for index in varobject:
                print (" ",index, varobject[index].value)

    def print_sol(self):
        for v in self.model.component_objects(pe.Var, active=True):

            # if str (v) =='var_time_a_d': 
            #     print ("Variable",v)
            #     varobject = getattr(self.model, str(v))
            #     for index in varobject:
            #         val  = varobject[index].value
            #         print (" ",index, val)

            if str (v) =='var_dlnk_path_occ': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    if val == 1.0:
                        print (" ",index, val,"  start,end dlnk %f %f, dv %f"%( 
                            self.model.par_t_start_dlnk[index[1],index[2]],
                            self.model.par_t_end_dlnk[index[1],index[2]],
                            self.model.par_dlnk_dv[index[1],index[2]])
                        )

            elif str (v) =='var_xlnk_path_occ': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    if val == 1.0:
                        print (" ",index, val,"  start,end xlnk %f %f, dv %f" %( 
                            self.model.par_t_start_xlnk[index[1],index[2],index[3]],
                            self.model.par_t_end_xlnk[index[1],index[2],index[3]],
                            self.model.par_xlnk_dv[index[1],index[2],index[3]])
                        ) 

            elif str (v) =='var_path_indic': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" ",index, val)

            elif str (v) =='var_path_indic_sat': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" ",index, val)
            
            elif str (v) =='var_path_dv': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" ",index, val)

            # elif str (v) =='var_route_fail': 
            #     print ("Variable",v)
            #     varobject = getattr(self.model, str(v))
            #     for index in varobject:
            #         val  = varobject[index].value
            #         print (" ",index, val)

            # elif str (v) =='var_path_dv_dlnk': 
            #     print ("Variable",v)
            #     varobject = getattr(self.model, str(v))
            #     for index in varobject:
            #         val  = varobject[index].value
            #         print (" ",index, val)

    def extract_routes( self,verbose = False):
        #  note that routes are the same as paths

        selected_routes_by_index = {}
        route_dv_by_index = {}
        senses =  {}

        try:
            # get the down links and index them by path
            for dlnk_subscr in self.model.dlnk_path_subscripts:
                # These subscript indices are defined above in get_downlink_model_objects
                p = dlnk_subscr[0]
                if self.model.var_dlnk_path_occ[dlnk_subscr] == 1.0:
                    if not p in selected_routes_by_index.keys ():
                        selected_routes_by_index[p] = []
                    sat_indx = dlnk_subscr[1]
                    dlnk_indx = dlnk_subscr[2]
                    dlnk_wind =  self.dlink_winds_flat[sat_indx][dlnk_indx]
                    senses[dlnk_wind] = dlnk_wind.sat_indx
                    selected_routes_by_index[p].append(dlnk_wind)

            #  get the cross-links and index them by path
            for xlnk_subscr in self.model.xlnk_path_subscripts:
                # These subscript indices are defined above in get_downlink_model_objects
                p = xlnk_subscr[0]
                if self.model.var_xlnk_path_occ[xlnk_subscr] == 1.0:
                    if not p in selected_routes_by_index.keys ():
                        selected_routes_by_index[p] = []
                    sat_indx = xlnk_subscr[1]
                    other_sat_indx = xlnk_subscr[2]
                    xlnk_indx = xlnk_subscr[3]

                    access_sat_indx = sat_indx
                    access_other_sat_indx = other_sat_indx
                    # remember that xlink_winds is upper triangular,  because it  would be symmetric across the diagonal. swap indexing to account for that
                    if sat_indx >other_sat_indx:
                        access_sat_indx=other_sat_indx
                        access_other_sat_indx = sat_indx

                    xlnk_wind =  self.xlink_winds[access_sat_indx][access_other_sat_indx][xlnk_indx]
                    #store the satellite from which data is being moved on this cross-link
                    senses[xlnk_wind] = sat_indx
                    selected_routes_by_index[p].append(xlnk_wind)

            #  add the observation
            for route in selected_routes_by_index. values ():
                senses[self.obs_wind] = self.obs_wind.sat_indx
                route.append (self.obs_wind)

            for p in self.model.paths:
                #  have to convert from pyomo  numeric type
                route_dv_by_index[p] =   pe.value(self.model.var_path_dv[p])

        except ValueError as e:
            print ("gp_route_selection.py: no routes could be extracted. problem is likely infeasible")

        #  now make the actual data route objects
        selected_routes =[]
        dr_id = 0
        for r in selected_routes_by_index. keys ():
            dr =  DataRoute (  
                ID=  dr_id, 
                route  = selected_routes_by_index[r],
                window_start_sats=senses,
                dv=route_dv_by_index[r]
            )
            dr.sort_windows ()
            selected_routes.append(dr)
            dr_id += 1

        if verbose:
            for dr in selected_routes:
                dr.print_route ( time_base = self.start_utc_dt) 

        return selected_routes