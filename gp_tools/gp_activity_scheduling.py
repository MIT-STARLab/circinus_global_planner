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

from .custom_activity_window import   ObsWindow,  DlnkWindow, XlnkWindow
from .routing_objects import DataRoute

class GPActivityScheduling():
    """docstring for GP activity scheduling"""
    def __init__(self,params):
        self.solver_max_runtime =params['solver_max_runtime_s']
        self.solver_name =params['solver_name']
        self.solver_run_remotely =params['solver_run_remotely']

    def get_stats(self,verbose=True):
        stats = {}
        stats['num_acts'] = sum([len ( self.all_acts_indcs)])
        stats['num_obs_acts'] = sum([len ( self.obs_act_indcs)])
        stats['num_link_acts'] = sum([len ( self.link_act_indcs)])
        stats['num_variables'] = self.model.nvariables ()
        stats['num_nobjectives'] = self.model.nobjectives ()
        stats['num_nconstraints'] = self.model.nconstraints ()

        if verbose:
            print ( "Considering %d activity windows" % (stats['num_acts']))
            print ( "Considering %d observation windows" % (stats['num_obs_acts']))
            print ( "Considering %d link windows" % (stats['num_link_acts']))
            print ( 'self.model.nvariables ()')
            print ( self.model.nvariables ())
            print ( 'self.model.nobjectives ()')
            print ( self.model.nobjectives ())
            print ( 'self.model.nconstraints ()')
            print ( self.model.nconstraints ())

        return stats

    def get_activity_structs( self,routes_flat):

        #  all activities are uniquely indexed. these structures keep track of those, and the mapping to activity objects
        all_acts_indcs = []
        #  these structures are for lookup in both directions
        all_acts_by_indx = {}
        all_acts_by_obj = {}

        # these structures keep track of the subset of unique indices that correspond to observations and links. Also we keep track of what data routes correspond to an activity
        link_act_indcs = []
        obs_act_indcs = []
        path_indcs_by_link_act = {}
        path_indcs_by_obs_act = {}
        dv_by_link_act = {}
        dv_by_obs_act = {}

        new_act_indx = 0
        for dr_indx, dr in enumerate (routes_flat):
            for act in dr.route:

                # if we haven't yet seen this activity, then add it to bookkeeping
                if not act in all_acts_by_obj.keys():
                    act_indx = new_act_indx
                    all_acts_indcs.append(act_indx)
                    all_acts_by_indx[act_indx] = act
                    all_acts_by_obj[act] = act_indx
                    new_act_indx += 1

                    # also need to add it to the list and dictionary for observations
                    if type(act) == ObsWindow:
                        obs_act_indcs.append(act_indx)
                        path_indcs_by_obs_act[act_indx] = []
                        path_indcs_by_obs_act[act_indx].append (dr_indx)
                        dv_by_obs_act[act_indx] = act.data_vol

                    # also need to add it to the list and dictionary for links
                    if type(act) == DlnkWindow or type(act) == XlnkWindow:
                        link_act_indcs.append(act_indx)
                        path_indcs_by_link_act[act_indx] = []
                        path_indcs_by_link_act[act_indx].append (dr_indx)
                        dv_by_link_act[act_indx] = act.data_vol

                #  if we have already seen the activity,  then just need to update the appropriate structures
                else:
                    act_indx = all_acts_by_obj[act]

                    # add the data route index
                    if type(act) == ObsWindow:
                        path_indcs_by_obs_act[act_indx].append (dr_indx)
                    if type(act) == DlnkWindow or type(act) == XlnkWindow:
                        path_indcs_by_link_act[act_indx].append (dr_indx)

        # print (path_indcs_by_link_act)
        # print (path_indcs_by_obs_act)

        return all_acts_indcs,obs_act_indcs,path_indcs_by_obs_act,dv_by_obs_act,link_act_indcs,path_indcs_by_link_act,dv_by_link_act
                    




    def make_model ( self,routes_flat, verbose = True):
        model = pe.ConcreteModel()

        # print(routes_flat)

        self.routes_flat = routes_flat

        if verbose:
            pass

        ##############################
        #  Make indices/ subscripts
        ##############################

        # note: a path is the same as a route

        (all_acts_indcs,
            obs_act_indcs,
            path_indcs_by_obs_act,
            dv_by_obs_act,
            link_act_indcs,
            path_indcs_by_link_act,
            dv_by_link_act) =  self.get_activity_structs(routes_flat)

        self.all_acts_indcs = all_acts_indcs
        self.obs_act_indcs = obs_act_indcs
        self.link_act_indcs = link_act_indcs

        #  subscript for each path p
        model.paths = pe.Set(initialize= range(len(routes_flat)))
        #  subscript for each activity a
        model.acts = pe.Set(initialize= all_acts_indcs)
        #  unique indices for observation and link acts
        model.obs_acts = pe.Set(initialize= obs_act_indcs)
        model.link_acts = pe.Set(initialize= link_act_indcs)

        ##############################
        #  Make parameters
        ##############################

        model.par_obs_dv = pe.Param(model.obs_acts,initialize =dv_by_obs_act)
        model.par_link_dv = pe.Param(model.link_acts,initialize =dv_by_link_act)
        model.par_path_dv = pe.Param(model.paths,initialize ={ i: dr.data_vol for i,dr in enumerate (routes_flat)})
        # each of these is essentially a dictionary indexed by link or obs act indx, with  the value being a list of path indices that are included within that act
        # these are valid indices into model.paths
        model.par_path_subscrs_by_link_act = pe.Param(model.link_acts,initialize =path_indcs_by_link_act)
        model.par_path_subscrs_by_obs_act = pe.Param(model.obs_acts,initialize =path_indcs_by_obs_act)

        ##############################
        #  Make variables
        ##############################


        # activity utilization variable indicating how much of an activity's capacity is used [1]
        model.var_activity_utilization  = pe.Var (model.acts, bounds =(0,1))
        # path utilization variable indicating how much of a path's capacity is used [2]
        model.var_path_utilization  = pe.Var (model.paths, bounds =(0,1))
        

        ##############################
        #  Make constraints
        ##############################

        # TODO: renumber  these with the final numbering

        def c1_rule( model,a):
            return (model.par_link_dv[a]*model.var_activity_utilization[a] -
                        sum(model.par_path_dv[p]*model.var_path_utilization[p] 
                            for p in model.par_path_subscrs_by_link_act[a]) 
                    >= 0)
        model.c1 =pe.Constraint ( model.link_acts,  rule=c1_rule)

        def c3b_rule( model,o):
            return (sum(model.par_path_dv[p]*model.var_path_utilization[p] 
                        for p in model.par_path_subscrs_by_obs_act[o]) 
                   <= model.par_obs_dv[o])
        model.c3b =pe.Constraint ( model.obs_acts,  rule=c3b_rule)

        

        ##############################
        #  Make objective
        ##############################

        def obj_rule(model):
            return (
                # model.par_obj_weight1 * 1/self.num_paths/model.par_obs_dv * sum(model.var_path_dv_dlnk[p,i,k] for i,k in model.dlnk_subscripts for p in model.paths )  +
                # model.par_obj_weight2 * 1/self.num_paths                  * sum(model.var_path_indic[p] for p in model.paths) +
                # model.par_obj_weight3 * 1/self.num_paths                  * sum(model.var_dlnk_path_occ[p,i,k]*model.par_dlnk_sf[i,k] for i,k in model.dlnk_subscripts for p in model.paths )
                sum(model.par_path_dv[p]*model.var_path_utilization[p] for p in model.paths)
            )
        model.obj = pe.Objective( rule=obj_rule, sense=pe.maximize )

        self.model = model

    # taken in part from Jeff Menezes' code at https://github.mit.edu/jmenezes/Satellite-MILP/blob/master/sat_milp_pyomo.py
    def solve(self):

        solver = po.SolverFactory(self.solver_name)
        solver.options['timelimit'] = self.solver_max_runtime

        # if we're running things remotely, then we will use the NEOS server (https://neos-server.org/neos/)
        if self.solver_run_remotely:
            solver_manager = po.SolverManagerFactory('neos')
            results = solver_manager.solve(self.model, opt= solver)
        else:
            # tee=True displays solver output in the terminal
            # keepfiles=True  keeps files passed to and from the solver
            results =  solver.solve(self.model, tee=True, keepfiles= False)

        if (results.solver.status == po.SolverStatus.ok) and (results.solver.termination_condition == po.TerminationCondition.optimal):
            print('this is feasible and optimal')
        elif results.solver.termination_condition == po.TerminationCondition.infeasible:
            print ('infeasible')
        else:
            # something else is wrong
            print (results.solver)

    def print_sol_all(self):
        for v in self.model.component_objects(pe.Var, active=True):
            print ("Variable",v)
            varobject = getattr(self.model, str(v))
            for index in varobject:
                print (" ",index, varobject[index].value)

    def print_sol(self):
        for v in self.model.component_objects(pe.Var, active=True):
            if str (v) =='var_activity_utilization': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" ",index, val)

            elif str (v) =='var_path_utilization': 
                print ("Variable",v)
                varobject = getattr(self.model, str(v))
                for index in varobject:
                    val  = varobject[index].value
                    print (" ",index, val)
            
