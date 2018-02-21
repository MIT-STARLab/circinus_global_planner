# Contains functionality for turning input data structures into the 
# objects used by the global planner  scheduling module
# 
# 
# @author Kit Kennedy
#

from pyomo import environ  as pe
from pyomo import opt  as po

class GPDataPathSelection():
    """docstring for GPInputProcessor"""
    def __init__(self,params):
        self.num_sats=params['num_sats']
        self.num_paths=params['path_selection_num_paths']
        self.start_utc_dt  =params['start_utc_dt']
        self.end_utc_dt  =params['end_utc_dt']
        self.M_t_s= 1000000 #  ~11.6 days
        self.M_dv_Mb= 1000000000 #  one petabit
        self.min_path_dv =params['path_selection_min_path_dv']

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


        #  now loop through and include paths as an index
        tuple_list_paths = []

        for path_indx in  range (num_paths):
            for sat_indx in  range (num_sats):

                    for dlnk_indx in  range (len (dlink_winds_flat[sat_indx])):
                        tuple_list_paths.append ((path_indx,sat_indx,dlnk_indx))

        return tuple_list, tuple_list_paths, t_start_dict_by_tuple, t_end_dict_by_tuple, dv_dict_by_tuple


    def make_model ( self,obs_wind,dlink_winds_flat,xlink_winds):
        model = pe.ConcreteModel()

        ##############################
        #  Make indices/ subscripts
        ##############################

        # get lists of indices ( as tuples) and start and end times for downlinks
        (dlnk_subscripts, 
            dlnk_path_subscripts,
            dlnk_t_start_dict, 
            dlnk_t_end_dict,
            dlnk_dv_dict)  = self.get_downlink_model_objects(
                                        dlink_winds_flat,
                                        self.num_paths,
                                        self.num_sats,
                                        self.start_utc_dt)

        # get lists of indices ( as tuples) and start and end times for  crosslinks
        (xlnk_subscripts, 
            xlnk_path_subscripts,
            xlnk_t_start_dict, 
            xlnk_t_end_dict,
            xlnk_dv_dict)  = self.get_crosslink_model_objects(
                                        xlink_winds,
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

        model.par_t_start_xlnk =pe.Param ( model.xlnk_subscripts,initialize =xlnk_t_start_dict)
        model.par_t_end_xlnk =pe.Param ( model.xlnk_subscripts,initialize =xlnk_t_end_dict)

        model.par_t_end_obs = pe.Param (initialize=(obs_wind.end - self.start_utc_dt).total_seconds ())
        model.par_obs_occ = pe.Param(model.sats,initialize ={ i: int (obs_wind.sat_indx == i) for i in  model.sats})

        model.par_dlnk_dv = pe.Param (model.dlnk_subscripts,initialize =dlnk_dv_dict)
        model.par_xlnk_dv = pe.Param (model.xlnk_subscripts,initialize =xlnk_dv_dict)

        model.par_obs_dv = pe.Param (initialize = obs_wind.data_vol)

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
        #  indicator that path p has been selected I_p_i [7]
        model.var_path_indic_sat  = pe.Var (model.paths,  model.sats, within =pe.Binary)
        #  arrival and departure times for each path p through sat i  t_(p,i,a), t_(p,i,d) [5]
        model.var_time_a_d   = pe.Var (model.paths,  model.sats,  model.arrive_depart, within =pe.NonNegativeReals)


        # print (type (model.var_xlnk_path_occ))
        # print (model.var_xlnk_path_occ[ (1,1,2,1)])
        # print (model.var_time_a_d[ (0,1,'a')])
        # print ( type (model.var_xlnk_path_occ[ (1,1,2,1)]))

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

        #  replaced by  constraint 19
        # model.c8a = pe.ConstraintList()
        # for p in model.paths: 
        #     for xlnk_subsc in model.xlnk_subscripts:
        #         #  order of elements in subscript tuple is  (see get_crosslink_model_objects()):
        #         #    i          j          k
        #         # (sat_indx, xsat_indx,xlnk_indx)
        #         i = xlnk_subsc[0]
        #         j = xlnk_subsc[1]
        #         k = xlnk_subsc[2]

        #         #  the constraint itself  is symmetric in i and j,  so don't want to consider when j is less than i
        #         if j > i:
        #             model.c8a.add( model.var_xlnk_path_occ[p,i,j,k] + model.var_xlnk_path_occ[p,j,i,k] <= 1 )
        
        def c8b_rule( model,i,j,k):
            return  sum(model.var_xlnk_path_occ[p,i,j,k] for p in  model.paths) <= 1
        model.c8b =pe.Constraint ( model.xlnk_subscripts,  rule=c8b_rule)

        #  this constraint is taken over by 9b.  also it's wrong as written. j also should be able to be less than i.  forget the nonsymmetric crap
        # #  for constraint 9 
        # def  get_sat_pair_tuples_non_sym(num_sats):
        #     tuple_list =[]
        #     for i in range( num_sats):
        #         for j in range(i +1, num_sats):
        #             tuple_list.append ( (i,j))
        #     return  tuple_list

        #  same as ^ ^ this constraint is taken over by 9b.  also it's wrong as written. j also should be able to be less than i.  forget the nonsymmetric crap
        # # constraint 9
        # sat_non_sym_tuples =  get_sat_pair_tuples_non_sym (self.num_sats)        
        # model.c9  = pe.ConstraintList()
        # for p in model.paths: 
        #     for tup in sat_non_sym_tuples:
        #         # tup is  (i,j)
        #         i =  tup[0]
        #         j =  tup[1]
        #         model.c9.add( 
        #             sum ([
        #                 model.var_xlnk_path_occ[p,i,j,k]  + 
        #                 model.var_xlnk_path_occ[p,j,i,k] 
        #                 for k in range (len (xlink_winds[i][j]))
        #             ])  
        #             <= 1
        #         )

        #  replaced by 19
        # model.c9b  = pe.ConstraintList()
        # for p in range(1):  #model.paths: 
        #     for j in model.sats:
        #         j_terms = []

        #         # have to iterate over i in i
        #         for i in model.sats:
        #             if i == j:
        #                 continue

        #             num_i_j_xlnks = max(len(xlink_winds[i][j]),len(xlink_winds[j][i]))
        #             # for k in range (num_i_j_xlnks):
        #             #     print ('oui')
        #             #     j_terms.append (model.var_xlnk_path_occ[p,i,j,k])
        #             #     j_terms.append (model.var_xlnk_path_occ[p,j,i,k])

        #             j_terms += [
        #                     model.var_xlnk_path_occ[p,i,j,k]  + 
        #                     model.var_xlnk_path_occ[p,j,i,k] 
        #                     for k in range (num_i_j_xlnks)
        #                 ]

        #             # print ( model.var_xlnk_path_occ[p,i,j,0])
        #             # print ( type (model.var_xlnk_path_occ[p,i,j,0]))
        #             # print (i,j)
        #             # print ( [model.var_xlnk_path_occ[p,i,j,0]])
        #             # print ( [model.var_xlnk_path_occ[p,i,j,0],model.var_xlnk_path_occ[p,i,j,0]])
        #             # print ( [] + [model.var_xlnk_path_occ[p,i,j,0],model.var_xlnk_path_occ[p,i,j,0]])
        #             # print ( sum ([model.var_xlnk_path_occ[p,i,j,0],model.var_xlnk_path_occ[p,i,j,0]]))
        #             # print (j_terms)
        #             # print (sum ([
        #             #         model.var_xlnk_path_occ[p,i,j,k]  + 
        #             #         model.var_xlnk_path_occ[p,j,i,k] 
        #             #         for k in range (len (xlink_winds[i][j]))
        #             #     ]))

        #         # print (j_terms)
        #         # print (sum (j_terms))
        #         # input ()
        #         model.c9b.add( sum(j_terms) <= 1)


        # constraint 10
        model.c10  = pe.ConstraintList()
        for p in model.paths: 
            for j in  model.sats:
                for k_d in range(len (dlink_winds_flat[j])):

                    k_d_terms = []
                    for i in model.sats:
                        if i == j:
                            continue
                        
                        k_d_terms += [
                            model.var_xlnk_path_occ[p,i,j,k_x] 
                            for k_x in range (len (xlink_winds[i][j]))
                        ] 

                    model.c10.add(model.var_dlnk_path_occ[p,j,k_d]  <= sum(k_d_terms))


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

        # this constraint is not useful
        # def c14_rule( model,p,i,k):
        #     # return model.var_path_dv[p] <=  model.var_dlnk_path_occ[tuple([p,i,k])]*model.par_dlnk_dv[tuple([i,k])]
        #     # return model.var_path_dv[p] <=  model.par_dlnk_dv[tuple([i,k])]  +  self.M_dv_Mb*(1-model.var_dlnk_path_occ[tuple([p,i,k])])
        #     return model.var_path_dv[p] <= 2000*model.var_dlnk_path_occ[tuple([p,i,k])]
        # model.c14 =pe.Constraint (   model.dlnk_path_subscripts,  rule =c14_rule )

        #  this constraint is also deprecated.  redundant with 15b  (though this may actually be more efficient than 15b)
        # def c15_rule( model,p):
        #     return model.var_path_dv[p]  >= model.par_min_path_dv*model.var_path_indic[p]
        #     # return model.var_path_dv[p]  >= 10*model.var_path_indic[p]
        # model.c15 =pe.Constraint ( model.paths, rule=c15_rule) 

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

        #  constraint 19
        model.c19  = pe.ConstraintList()
        for p in model.paths: 
            for j in model.sats:
                j_terms = []

                #  could have been a list comprehension but whatever
                for k_d in range(len (dlink_winds_flat[j])):
                    j_terms.append (model.var_dlnk_path_occ[p,j,k_d] )

                for i in model.sats:
                    if i == j:
                        continue

                    #  use this to account for the fact that the cross-link Windows matrix is upper triangular
                    num_i_j_xlnks = max(len(xlink_winds[i][j]),len(xlink_winds[j][i]))

                    j_terms += [
                        model.var_xlnk_path_occ[p,i,j,k]  + 
                        model.var_xlnk_path_occ[p,j,i,k] 
                        for k in range (num_i_j_xlnks)
                    ]

                j_terms.append (model.par_obs_occ[j])

                model.c19.add( sum(j_terms) == 2*model.var_path_indic_sat[p,j])

        ##############################
        #  Make objective
        ##############################

        def obj_rule(model):
            # return sum( model.var_path_indic[p] for p in model.paths)
            return sum (model.var_path_dv_dlnk[p,i,k] for i,k in model.dlnk_subscripts for p in model.paths )
            # return sum( model.var_dlnk_path_occ[dlnk_subsc] for dlnk_subsc in  model.dlnk_path_subscripts)
        model.obj = pe.Objective( rule=obj_rule, sense=pe.maximize )

 
        
        
        
        

        self.model = model

    # lifted from  Jeff Menezes' code at https://github.mit.edu/jmenezes/Satellite-MILP/blob/master/sat_milp_pyomo.py
    def solve(self):
        solver = po.SolverFactory('gurobi')
        results =  solver.solve(self.model, tee=True, keepfiles=False, options_string=" time_limit=10")
        
        if (results.solver.status == po.SolverStatus.ok) and (results.solver.termination_condition == po.TerminationCondition.optimal):
            print('this is feasible and optimal')
        elif results.solver.termination_condition == po.TerminationCondition.infeasible:
            print ('infeasible')
        else:
            # something else is wrong
            print (results.solver)


        print ( 'self.model.nvariables ()')
        print ( self.model.nvariables ())
        print ( 'self.model.nobjectives ()')
        print ( self.model.nobjectives ())
        print ( 'self.model.nconstraints ()')
        print ( self.model.nconstraints ())
        # elif (results.solver.status != po.SolverStatus.ok):
        #     print('Check solver not ok?')
        # elif (results.solver.termination_condition != po.TerminationCondition.optimal):  
        #     print('Check solver optimality?') 

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
                        print (" ",index, val,"  start dlnk %f, dv %f"%( 
                            self.model.par_t_start_dlnk[index[1],index[2]],
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

            # elif str (v) =='var_path_dv_dlnk': 
            #     print ("Variable",v)
            #     varobject = getattr(self.model, str(v))
            #     for index in varobject:
            #         val  = varobject[index].value
            #         print (" ",index, val)
