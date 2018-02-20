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

    @staticmethod
    def get_crosslink_model_objects(xlink_winds,num_paths,num_sats,start_utc_dt):
        tuple_list =[]
        #  start and end times of cross-links, indexed by xlnk i,j,k tuple
        t_start_dict_by_tuple =  {}
        t_end_dict_by_tuple = {}

        # loop through all satellite index combinations i,j. Note that there is duplication here; i,j and j,i will refer to the same window, and produce different tuples. This is a feature not a bug. We need both directions to have both directions represented
        # explicitly indexing by num_sats just to include a bit of error checking
        for sat_indx in  range (num_sats):
            for xsat_indx in  range (num_sats):
                if xsat_indx == sat_indx:
                    continue

                for xlnk_indx in  range (len (xlink_winds[sat_indx][xsat_indx])):
                    tup = (sat_indx, xsat_indx,xlnk_indx)
                    tuple_list.append (tup)
                    wind = xlink_winds[sat_indx][xsat_indx][xlnk_indx]
                    t_start_dict_by_tuple[tup] =  (wind.start - start_utc_dt).total_seconds ()
                    t_end_dict_by_tuple[tup] = (wind.end - start_utc_dt).total_seconds ()


        #  now loop through and include paths as an index
        tuple_list_paths = []

        for path_indx in  range (num_paths):
            for sat_indx in  range (num_sats):
                for xsat_indx in  range (num_sats):
                    if xsat_indx == sat_indx:
                        continue

                    for xlnk_indx in  range (len (xlink_winds[sat_indx][xsat_indx])):
                        tuple_list_paths.append ((path_indx,sat_indx, xsat_indx,xlnk_indx))

        return tuple_list, tuple_list_paths, t_start_dict_by_tuple, t_end_dict_by_tuple


    @staticmethod
    def  get_dlnk_path_subscripts(num_paths,num_sats,sats_num_dlnks=[]):
        tuple_list = []

        for path_indx in  range (num_paths):
            for sat_indx in  range (num_sats):
                    for dlnk_indx in  range (sats_num_dlnks[sat_indx]):
                        tuple_list.append ((path_indx,sat_indx,dlnk_indx))

        return tuple_list


    @staticmethod
    def  get_dlnk_subscripts(num_sats,sats_num_dlnks=[]):
        tuple_list = []

        for sat_indx in  range (num_sats):
            for dlnk_indx in  range (sats_num_dlnks[sat_indx]):
                tuple_list.append  ((sat_indx,dlnk_indx))

        return tuple_list


    def make_model ( self,obs,dlink_winds,xlink_winds):
        model = pe.ConcreteModel()

        ##############################
        #  Make indices/ subscripts
        ##############################

        # create lists of tuples that serve as subscripts for cross-links and down links
        # sats_num_xlnks = [len(sat_winds) for sat_winds in xlink_winds]
        sats_num_dlnks = [len(sat_winds) for sat_winds in dlink_winds]
        # xlnk_path_subscripts = self.get_crosslink_path_subscripts( self.num_paths, self.num_sats,sats_num_xlnks)
        dlnk_path_subscripts = self.get_dlnk_path_subscripts( self.num_paths, self.num_sats,sats_num_dlnks)
        dlnk_subscripts = self.get_dlnk_subscripts( self.num_sats,sats_num_dlnks)
        # xlnk_subscripts = self.get_xlnk_subscripts( self.num_sats,sats_num_xlnks)

        (xlnk_subscripts, 
            xlnk_path_subscripts,
            xlnk_t_start_dict, 
            xlnk_t_end_dict)  = self.get_crosslink_model_objects(
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

        # model.pprint ()
        print (len(xlnk_path_subscripts))
        print (len(dlnk_path_subscripts))
        print (len(dlnk_subscripts))

        ##############################
        #  Make parameters
        ##############################

        #  start from 0  for our time system
        model.par_t_start = pe.Param (initialize=0)
        model.par_t_end = pe.Param (initialize= (self.end_utc_dt - self.start_utc_dt).total_seconds ())

        model.par_t_end_xlnk =pe.Param ( model.xlnk_subscripts,initialize =xlnk_t_end_dict)

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
        #  indicator that path p is been selected I_p [4]
        model.var_path_indic  = pe.Var (model.paths, within =pe.Binary)
        #  arrival and departure times for each path p through sat i  t_(p,i,a), t_(p,i,d) [5]
        model.var_time_a_d   = pe.Var (model.paths,  model.sats,  model.arrive_depart, within =pe.NonNegativeReals)


        # print (type (model.var_xlnk_path_occ))
        # print (model.var_xlnk_path_occ[ (1,1,2,1)])
        # print (model.var_time_a_d[ (0,1,'a')])
        # print ( type (model.var_xlnk_path_occ[ (1,1,2,1)]))

        ##############################
        #  Make constraints
        ##############################

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

        # def c3_rule( model,p,i):
        #      return model.var_time_a_d [p,i,'a']  <= model.var_time_a_d [p,i,'d']
        # model.c3 =pe.Constraint ( model.paths, model.sats,  rule=c3_rule)

 



        ##############################
        #  Make objective
        ##############################

        def obj_rule(model):
            return sum( model.var_time_a_d [p,i,'a'] for p in model.paths for i in model.sats)
        model.obj = pe.Objective( rule=obj_rule, sense=pe.maximize )

        self.model = model

    # lifted from  Jeff Menezes' code at https://github.mit.edu/jmenezes/Satellite-MILP/blob/master/sat_milp_pyomo.py
    def solve(self):
        solver = po.SolverFactory('gurobi')
        results =  solver.solve(self.model, tee=True, keepfiles=False, options_string="mip_tolerances_integrality=1e-9 mip_tolerances_mipgap=0")
        
        if (results.solver.status == po.SolverStatus.ok) and (results.solver.termination_condition == po.TerminationCondition.optimal):
            print('this is feasible and optimal')
        elif results.solver.termination_condition == po.TerminationCondition.infeasible:
            print ('infeasible')
        else:
            # something else is wrong
            print (results.solver)
        # elif (results.solver.status != po.SolverStatus.ok):
        #     print('Check solver not ok?')
        # elif (results.solver.termination_condition != po.TerminationCondition.optimal):  
        #     print('Check solver optimality?') 

    # lifted from  Jeff Menezes' code at https://github.mit.edu/jmenezes/Satellite-MILP/blob/master/sat_milp_pyomo.py
    def print_sol(self):
        for v in self.model.component_objects(pe.Var, active=True):
            print ("Variable",v)
            varobject = getattr(self.model, str(v))
            for index in varobject:
                print (" ",index, varobject[index].value)
