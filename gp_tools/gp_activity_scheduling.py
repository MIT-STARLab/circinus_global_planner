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

from .routing_objects import DataRoute

class GPActivityScheduling():
    """docstring for GP route selection"""
    def __init__(self,params):
        pass

    def make_model ( self,obs_wind,dlink_winds_flat,xlink_winds, verbose = True):
        model = pe.ConcreteModel()

        if verbose:
            pass

        ##############################
        #  Make indices/ subscripts
        ##############################

        # ge
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


        ##############################
        #  Make variables
        ##############################

        #  path occupancy variables  for each crosslink x_(p,i,j,k) [1]
        model.var_xlnk_path_occ  = pe.Var (model.xlnk_path_subscripts, within =pe.Binary)
        


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

    # lifted from  Jeff Menezes' code at https://github.mit.edu/jmenezes/Satellite-MILP/blob/master/sat_milp_pyomo.py
    def solve(self):
        solver = po.SolverFactory('gurobi')
        results =  solver.solve(self.model, tee=True, keepfiles=False, options_string=" time_limit=%f" % ( self.solver_max_runtime))
        
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