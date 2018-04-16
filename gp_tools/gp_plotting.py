#  plotting tools for global planner
# 
# @author Kit Kennedy

import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig
from matplotlib.patches import Rectangle, Circle
import numpy as np

from .routing_objects import DataRoute

class GPPlotting():

    xlnk_colors = [
        '#FF0000',
        '#FF3399',
        '#990000',
        '#990099',
        '#FF9900'
    ]
    xlnk_color_rollover = len(xlnk_colors)

    # hard coding that we will use the first data route index in the list of data indices to decide the color for a link
    route_index_to_use = 0

    """docstring for GP plotting"""
    def __init__(self,gp_params):
        """initializes based on parameters
        
        initializes based on parameters
        :param gp_params: global namespace parameters created from input files (possibly with some small non-structural modifications to params). The name spaces here should trace up all the way to the input files.
        :type params: dict
        """
        plot_params = gp_params['gp_general_params']['plot_params']
        sat_params = gp_params['gp_orbit_prop_params']['sat_params']
        self.obs_params = gp_params['gp_orbit_prop_params']['obs_params']

        self.plot_fig_extension=plot_params['plot_fig_extension']
        self.time_units=plot_params['time_units']
        self.winds_plot_obs=plot_params['winds_plot_obs']
        self.winds_plot_obs_choices=plot_params['winds_plot_obs_choices']
        self.winds_plot_dlnks=plot_params['winds_plot_dlnks']
        self.winds_plot_dlnks_choices=plot_params['winds_plot_dlnks_choices']
        self.winds_plot_xlnks=plot_params['winds_plot_xlnks']
        self.winds_plot_xlnks_choices=plot_params['winds_plot_xlnks_choices']
        self.energy_usage_plot_params=plot_params['energy_usage_plot_params']
        self.data_usage_plot_params=plot_params['data_usage_plot_params']
        self.obs_aoi_metrics_plot_params=plot_params['obs_aoi_metrics_plot_params']
        self.cmd_aoi_metrics_plot_params=plot_params['cmd_aoi_metrics_plot_params']
        self.tlm_aoi_metrics_plot_params=plot_params['tlm_aoi_metrics_plot_params']

        self.sat_id_order = sat_params['sat_id_order']
        self.all_targ_IDs = [targ['id'] for targ in self.obs_params['targets']]

    def plot_winds(
        self,
        sats_ids_list,
        all_obs_winds_flat,
        obs_winds_flat,
        all_dlnk_winds_flat,
        dlnk_winds_flat, 
        all_xlnk_winds_flat,
        xlnk_winds_flat,
        route_indcs_by_wind,
        plot_start,
        plot_end,
        plot_title = 'Route Plot', 
        plot_size_inches = (12,12),
        plot_include_dlnk_labels = True,
        plot_include_xlnk_labels = True,
        show=False,
        fig_name='plots/xlnk_dlnk_plot.pdf'):
        '''
        Displays a 2D plot of assignments for each agent with respect to time

        '''

        if self.time_units == 'hours':
            time_divisor = 3600
        if self.time_units == 'minutes':
            time_divisor = 60
        
        time_to_end = (plot_end-plot_start).total_seconds()/time_divisor

        num_sats = len(sats_ids_list)

        #  make a new figure
        plt.figure()

        #  create subplots for satellites
        fig = plt.gcf()
        fig.set_size_inches( plot_size_inches)
        # print fig.get_size_inches()

        plt.title( plot_title)

        # have to sort before applying axis labels, otherwise x label shows up in a weird place
        # sats.sort(key=lambda x: x.agent_ID)

        # keep a running list of all the window IDs seen,  which we'll use for a sanity check
        all_wind_ids = []

        # for each agent
        for  plot_indx, sat_id in enumerate (sats_ids_list):
            #  get the index for this ID
            sat_indx = self.sat_id_order.index(str(sat_id))

            # 
            axes = plt.subplot( num_sats,1,plot_indx+1)
            if plot_indx == np.floor(num_sats/2):
                plt.ylabel('Satellite Index\n\n' + str(sat_indx))
            else:
                plt.ylabel('' + str(sat_indx))


            # no y-axis labels
            plt.tick_params(
                axis='y',
                which='both',
                left='off',
                right='off',
                labelleft='off'
            )

            SMALL_SIZE = 8
            MEDIUM_SIZE = 10
            BIGGER_SIZE = 15

            plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
            plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
            plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
            plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
            plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
            plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

            # set axis length. Time starts at 0
            plt.axis((0, time_to_end, 0, 2))

            current_axis = plt.gca()

            #  this is used to alternate the vertical position of the cross-link rectangle
            obs_rectangle_rotator = 0
            obs_choices_rectangle_rotator = 0
            obs_rotation_rollover = 2
            #  this is used to alternate the vertical position of the cross-link rectangle
            xlnk_rectangle_rotator = 0
            xlnk_choices_rectangle_rotator = 0
            xlnk_rotation_rollover = 2
            #  this is used to alternate the vertical position of labels
            xlnk_label_rotator = 0
            
            #  this is used to alternate the vertical position of the  downlink rectangle
            dlnk_rectangle_rotator = 0
            dlnk_choices_rectangle_rotator = 0
            dlnk_rotation_rollover = 2
            dlnk_label_rotator = 0

            #  these hold the very last plot object of a given type added. Used for legend below
            d_w = None
            d = None
            x_w = None
            x = None

            ##########################
            # plot window choices
            ##########################

            obs_choices_rectangle_rotator_hist = {}
            dlnk_choices_rectangle_rotator_hist = {}
            xlnk_choices_rectangle_rotator_hist = {}

            # plot the crosslink "choices" -  meant to represent the windows that could have been chosen
            #  plot cross-links first, so that they are the furthest back (lowest z value) on the plot, and observations and downlinks will appear on top ( because there are generally a lot more cross-links than observations and down links)
            if self.winds_plot_xlnks_choices:
                num_xlnk_wind = 0

                if len(all_xlnk_winds_flat) > 0:
                    for xlnk_wind in all_xlnk_winds_flat[sat_indx]:

                        xlnk_wind_start = (xlnk_wind.start-plot_start).total_seconds()/time_divisor
                        xlnk_wind_end = (xlnk_wind.end-plot_start).total_seconds()/time_divisor

                        # plot the task duration
                        bottom_vert_loc = xlnk_choices_rectangle_rotator
                        x_w = Rectangle((xlnk_wind_start, bottom_vert_loc), xlnk_wind_end-xlnk_wind_start, bottom_vert_loc+1,alpha=1,fill=True,color='#FFBCBC')

                        current_axis.add_patch(x_w)

                        # plt.text(xlnk_wind_start+0.15, bottom_vert_loc+0.1, xlnk_wind.window_ID , fontsize=10, color = 'k')

                        # save off the rotator choice so that we can look it up again
                        xlnk_choices_rectangle_rotator_hist[xlnk_wind] = xlnk_choices_rectangle_rotator

                        xlnk_choices_rectangle_rotator =  (xlnk_choices_rectangle_rotator+1)%xlnk_rotation_rollover
                        num_xlnk_wind += 1

            # plot the observation "choices" -  meant to represent the windows that could have been chosen
            if self.winds_plot_obs_choices:
                if len(all_obs_winds_flat) > 0:
                    for obs_wind in all_obs_winds_flat[sat_indx]:

                        obs_start = (obs_wind.start-plot_start).total_seconds()/time_divisor
                        obs_end = (obs_wind.end-plot_start).total_seconds()/time_divisor

                        # plot the task duration
                        bottom_vert_loc = obs_choices_rectangle_rotator
                        d = Rectangle((obs_start, bottom_vert_loc), obs_end-obs_start, bottom_vert_loc+1,alpha=1,fill=True,color='#BFFFBF')
                        current_axis.add_patch(d)

                        plt.text(obs_end+0.15, bottom_vert_loc+0.1, obs_wind.target_IDs , fontsize=10, color = 'k')

                        # save off the rotator choice so that we can look it up again
                        obs_choices_rectangle_rotator_hist[obs_wind] = obs_choices_rectangle_rotator

                        obs_choices_rectangle_rotator =  (obs_choices_rectangle_rotator+1)%obs_rotation_rollover

            # plot the downlink "choices" -  meant to represent the windows that could have been chosen
            if self.winds_plot_dlnks_choices:
                num_dlnk_wind = 0

                if len(all_dlnk_winds_flat) > 0:
                    for dlnk_wind in all_dlnk_winds_flat[sat_indx]:

                        dlnk_wind_start = (dlnk_wind.start-plot_start).total_seconds()/time_divisor
                        dlnk_wind_end = (dlnk_wind.end-plot_start).total_seconds()/time_divisor

                        # plot the task duration
                        bottom_vert_loc = dlnk_choices_rectangle_rotator
                        d_w = Rectangle((dlnk_wind_start, bottom_vert_loc), dlnk_wind_end-dlnk_wind_start, bottom_vert_loc+1,alpha=1,fill=True,color='#BFBFFF')

                        current_axis.add_patch(d_w)
                        # plt.text( (dlnk_wind_end+dlnk_wind_start)/2 - 0.15, 0.1, dlnk_wind.gs_indx , fontsize=10, color = 'k')

                        # plt.text(dlnk_wind_start+0.15, bottom_vert_loc+0.1, dlnk_wind.window_ID , fontsize=10, color = 'k')

                        # save off the rotator choice so that we can look it up again
                        dlnk_choices_rectangle_rotator_hist[dlnk_wind] = dlnk_choices_rectangle_rotator

                        dlnk_choices_rectangle_rotator =  (dlnk_choices_rectangle_rotator+1)%dlnk_rotation_rollover

                        num_dlnk_wind += 1

            ##########################
            # plot windows executed
            ##########################

            #  plot the executed cross-links
            #  plot cross-links first, so that they are the furthest back (lowest z value) on the plot, and observations and downlinks will appear on top ( because there are generally a lot more cross-links than observations and down links)
            if self.winds_plot_xlnks:
                num_xlnk_exe = 0
                if len(xlnk_winds_flat) > 0:
                    for xlnk_wind in xlnk_winds_flat[sat_indx]:

                        xlnk_start = (xlnk_wind.start-plot_start).total_seconds()/time_divisor
                        xlnk_end = (xlnk_wind.end-plot_start).total_seconds()/time_divisor

                        #  update the rotator value if we've already added this window to the plot in the "choices" code above
                        if xlnk_wind in xlnk_choices_rectangle_rotator_hist.keys ():
                            xlnk_rectangle_rotator = xlnk_choices_rectangle_rotator_hist[xlnk_wind]

                        bottom_vert_loc = xlnk_rectangle_rotator

                        # if we have been given route indices, then then use one of them to configure the crosslink box color
                        dr_indx = None
                        if (len(route_indcs_by_wind.keys())) > 0:
                            dr_indx = route_indcs_by_wind[xlnk_wind][self.route_index_to_use]
                            xlnk_color_indx = dr_indx %  self.xlnk_color_rollover
                        #  otherwise just go with the first color
                        else:
                            xlnk_color_indx = 0

                        xlnk_color = self.xlnk_colors[xlnk_color_indx]
                        
                        # plot the task duration
                        x = Rectangle((xlnk_start, bottom_vert_loc), xlnk_end-xlnk_start, bottom_vert_loc +1,alpha=1,fill=False,color=xlnk_color,hatch='///////')
                        current_axis.add_patch(x)

                        xlnk_rectangle_rotator =  (xlnk_rectangle_rotator+1)%xlnk_rotation_rollover

                        if plot_include_xlnk_labels:
                            other_sat_indx = xlnk_wind.xsat_indx if xlnk_wind.xsat_indx != sat_indx else xlnk_wind.sat_indx

                            #   put label in desired vertical spot
                            left_horizontal_loc = xlnk_start + 0.15

                            #  again, if we know route indices include them in label
                            if dr_indx:
                                label_text = "%d,%d" %(dr_indx,other_sat_indx)
                            else:                                
                                label_text = "%d" %(other_sat_indx)

                            # update the rotator value if we've already added this window to the plot in the "choices" code above
                            if xlnk_wind in xlnk_choices_rectangle_rotator_hist.keys ():
                                xlnk_label_rotator = xlnk_choices_rectangle_rotator_hist[xlnk_wind]

                            if xlnk_label_rotator == 0:
                                plt.text(left_horizontal_loc, 0.1, label_text , fontsize=10, color = 'k')
                            elif xlnk_label_rotator == 1:
                                plt.text( left_horizontal_loc, 1.1, label_text , fontsize=10, color = 'k')

                            xlnk_label_rotator = (xlnk_label_rotator+1)%xlnk_rotation_rollover

            # plot the observations that are actually executed
            if self.winds_plot_obs:
                if len(obs_winds_flat) > 0:
                    for obs_wind in obs_winds_flat[sat_indx]:

                        obs_start = (obs_wind.start-plot_start).total_seconds()/time_divisor
                        obs_end = (obs_wind.end-plot_start).total_seconds()/time_divisor

                        #  update the rotator value if we've already added this window to the plot in the "choices" code above
                        if obs_wind in obs_choices_rectangle_rotator_hist.keys ():
                            obs_rectangle_rotator = obs_choices_rectangle_rotator_hist[obs_wind]

                        # plot the task duration
                        bottom_vert_loc = obs_rectangle_rotator
                        d = Rectangle((obs_start, bottom_vert_loc), obs_end-obs_start, bottom_vert_loc+1,alpha=1,fill=False,color='#00FF00',hatch='///////')
                        current_axis.add_patch(d)

                        obs_rectangle_rotator =  (obs_rectangle_rotator+1)%obs_rotation_rollover

                        if obs_wind.window_ID in all_wind_ids:
                            raise Exception('Found a duplicate unique window ID where it should not have been possible')
                        all_wind_ids.append(obs_wind.window_ID)

            # plot the executed down links
            if self.winds_plot_dlnks:
                num_dlnk_exe = 0
                if len(dlnk_winds_flat) > 0:
                    for dlnk_wind in dlnk_winds_flat[sat_indx]:

                        dlnk_start = (dlnk_wind.start-plot_start).total_seconds()/time_divisor
                        dlnk_end = (dlnk_wind.end-plot_start).total_seconds()/time_divisor

                        gs_indx = dlnk_wind.gs_indx

                        #  update the rotator value if we've already added this window to the plot in the "choices" code above
                        if dlnk_wind in dlnk_choices_rectangle_rotator_hist.keys ():
                            dlnk_rectangle_rotator = dlnk_choices_rectangle_rotator_hist[dlnk_wind]

                        # plot the task duration
                        bottom_vert_loc = dlnk_rectangle_rotator
                        d = Rectangle((dlnk_start, bottom_vert_loc), dlnk_end-dlnk_start, bottom_vert_loc+1,alpha=1,fill=False,color='#0000FF',hatch='///////')
                        current_axis.add_patch(d)

                        dlnk_rectangle_rotator =  (dlnk_rectangle_rotator+1)%dlnk_rotation_rollover

                        if dlnk_wind.window_ID in all_wind_ids:
                            raise Exception('Found a duplicate unique window ID where it should not have been possible')
                        all_wind_ids.append(dlnk_wind.window_ID)

                        if plot_include_dlnk_labels:
                            label_text = ""

                            # if we have been given route indices, then include them in the label
                            if (len(route_indcs_by_wind.keys())) > 0:
                                dr_indcs = route_indcs_by_wind[dlnk_wind]
                                for dr_indx in dr_indcs:
                                    label_text += "%d,"%(dr_indx)
                                label_text += ";"
                            
                            #  add the ground station index to the label
                            label_text += "%d"%(gs_indx)

                            #   put label in desired vertical spot
                            left_horizontal_loc = dlnk_start + 0.15
                            
                            #  update the rotator value if we've already added this window to the plot in the "choices" code above
                            if dlnk_wind in dlnk_choices_rectangle_rotator_hist.keys ():
                                dlnk_label_rotator = dlnk_choices_rectangle_rotator_hist[dlnk_wind]
                            
                            if dlnk_label_rotator == 0:
                                plt.text(left_horizontal_loc, 0.1, label_text , fontsize=10, color = 'k')
                            elif dlnk_label_rotator == 1:
                                plt.text( left_horizontal_loc, 1.1, label_text , fontsize=10, color = 'k')

                            dlnk_label_rotator = (dlnk_label_rotator+1)%dlnk_rotation_rollover

                        num_dlnk_exe += 1

            #  if were at the last satellite ( at the bottom of all the plots), then add X axis labels
            if not plot_indx+1 == num_sats:
                ax = plt.gca()
                plt.setp(ax.get_xticklabels(), visible=False)

        legend_objects = []
        legend_objects_labels = []
        if d_w: 
            legend_objects.append(d_w)
            legend_objects_labels.append('D all')
        if d: 
            legend_objects.append(d)
            legend_objects_labels.append('Dlnk')
        if x_w: 
            legend_objects.append(x_w)
            legend_objects_labels.append('X all')
        if x: 
            legend_objects.append(x)
            legend_objects_labels.append('Xlnk')

        plt.legend(legend_objects, legend_objects_labels ,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.xlabel('Time (%s)'%(self.time_units))

        # use the last axes to set the entire plot background color
        axes.patch.set_facecolor('w')

        if show:
            plt.show()
        else:
            savefig(fig_name,format=self.plot_fig_extension)

    def plot_data_circles(
        self,
        sats_ids_list,
        all_obs_winds_flat,
        obs_winds_flat,
        all_dlnk_winds_flat,
        dlnk_winds_flat,
        all_xlnk_winds_flat,
        xlnk_winds_flat,
        route_indcs_by_wind,
        plot_start,
        plot_end,
        plot_title = 'Route Data Plot', 
        plot_size_inches = (12,12),
        plot_include_labels = True,
        show=False,
        fig_name='plots/xlnk_dlnk_plot.pdf'):
        '''
        Displays a 2D plot of assignments for each agent with respect to time

        '''

        if self.time_units == 'hours':
            time_divisor = 3600
        if self.time_units == 'minutes':
            time_divisor = 60
        
        time_to_end = (plot_end-plot_start).total_seconds()/time_divisor

        num_sats = len(sats_ids_list)

        def dv_radius_scale(data_vol):
            scaler = 100*10
            radius = data_vol/scaler
            if radius > 1.5:
                radius = 1
            if radius < 0.1:
                radius = 0.1

            return radius

        #  make a new figure
        plt.figure()

        #  create subplots for satellites
        fig = plt.gcf()
        fig.set_size_inches( plot_size_inches)
        # print fig.get_size_inches()

        plt.title( plot_title)

        # have to sort before applying axis labels, otherwise x label shows up in a weird place
        # sats.sort(key=lambda x: x.agent_ID)

        # keep a running list of all the window IDs seen,  which we'll use for a sanity check
        all_wind_ids = []

        # for each agent
        for  plot_indx, sat_id in enumerate (sats_ids_list):
            #  get the index for this ID
            sat_indx = self.sat_id_order.index(str(sat_id))

            # 
            axes = plt.subplot( num_sats,1,plot_indx+1)
            if plot_indx == np.floor(num_sats/2):
                plt.ylabel('Satellite Index\n\n' + str(sat_indx))
            else:
                plt.ylabel('' + str(sat_indx))


            # no y-axis labels
            plt.tick_params(
                axis='y',
                which='both',
                left='off',
                right='off',
                labelleft='off'
            )

            # set axis length. Time starts at 0
            plt.axis((0, time_to_end, 0, 2))

            current_axis = plt.gca()

            #  this is used to alternate the vertical position of the cross-link rectangle
            obs_rectangle_rotator = 0
            obs_choices_rectangle_rotator = 0
            obs_rotation_rollover = 4
            #  this is used to alternate the vertical position of the cross-link rectangle
            xlnk_rectangle_rotator = 0
            xlnk_choices_rectangle_rotator = 0
            xlnk_rotation_rollover = 4
            #  this is used to alternate the vertical position of labels
            xlnk_label_rotator = 0
            
            #  this is used to alternate the vertical position of the  downlink rectangle
            dlnk_rectangle_rotator = 0
            dlnk_choices_rectangle_rotator = 0
            dlnk_rotation_rollover = 4
            dlnk_label_rotator = 0

            #  these hold the very last plot object of a given type added. Used for legend below
            d_w = None
            d = None
            x_w = None
            x = None

            ###################
            # obs
            ###################

            obs_choices_rectangle_rotator_hist = {}
            dlnk_choices_rectangle_rotator_hist = {}
            xlnk_choices_rectangle_rotator_hist = {}

            # plot the  observations
            if self.winds_plot_obs_choices:
                if len(all_obs_winds_flat) > 0:
                    for obs_wind in all_obs_winds_flat[sat_indx]:

                        obs_start = (obs_wind.start-plot_start).total_seconds()/time_divisor
                        obs_end = (obs_wind.end-plot_start).total_seconds()/time_divisor
                        middle = (obs_end + obs_start)/2

                        # plot the task duration
                        bottom_vert_loc = obs_choices_rectangle_rotator
                        radius =  dv_radius_scale(obs_wind.data_vol)
                        o_w = Circle((middle, bottom_vert_loc), radius,alpha=1,fill=True,color='#BFFFBF')
                        current_axis.add_patch(o_w)

                        # plt.text(obs_start+0.15, bottom_vert_loc+0.1, obs_wind.window_ID , fontsize=10, color = 'k')

                        # save off the rotator choice so that we can look it up again
                        obs_choices_rectangle_rotator_hist[obs_wind] = obs_choices_rectangle_rotator

                        obs_choices_rectangle_rotator =  (obs_choices_rectangle_rotator+2)%obs_rotation_rollover

            if self.winds_plot_obs:
                if len(obs_winds_flat) > 0:
                    for obs_wind in obs_winds_flat[sat_indx]:

                        obs_start = (obs_wind.start-plot_start).total_seconds()/time_divisor
                        obs_end = (obs_wind.end-plot_start).total_seconds()/time_divisor
                        middle = (obs_end + obs_start)/2

                        #  update the rotator value if we've already added this window to the plot in the "choices" code above
                        if obs_wind in obs_choices_rectangle_rotator_hist.keys ():
                            obs_rectangle_rotator = obs_choices_rectangle_rotator_hist[obs_wind]

                        # plot the task duration
                        bottom_vert_loc = obs_rectangle_rotator
                        radius =  dv_radius_scale(obs_wind.scheduled_data_vol)
                        o = Circle((middle, bottom_vert_loc), radius,alpha=1,fill=False,color='#00FF00',hatch='///////')
                        current_axis.add_patch(o)

                        obs_rectangle_rotator =  (obs_rectangle_rotator+2)%obs_rotation_rollover

                        if obs_wind.window_ID in all_wind_ids:
                            raise Exception('Found a duplicate unique window ID where it should not have been possible')
                        all_wind_ids.append(obs_wind.window_ID)

            ###################
            # dlnks
            ###################

            #  plot the potential down links
            if self.winds_plot_dlnks_choices:
                num_dlnk_wind = 0

                if len(all_dlnk_winds_flat) > 0:
                    for dlnk_wind in all_dlnk_winds_flat[sat_indx]:

                        dlnk_wind_start = (dlnk_wind.start-plot_start).total_seconds()/time_divisor
                        dlnk_wind_end = (dlnk_wind.end-plot_start).total_seconds()/time_divisor
                        middle = (dlnk_wind_end + dlnk_wind_start)/2

                        # plot the task duration
                        bottom_vert_loc = dlnk_choices_rectangle_rotator
                        radius =  dv_radius_scale(dlnk_wind.data_vol)
                        d_w = Circle((middle, bottom_vert_loc), radius,alpha=1,fill=True,color='#BFBFFF')

                        current_axis.add_patch(d_w)
                        # plt.text( (dlnk_wind_end+dlnk_wind_start)/2 - 0.15, 0.1, dlnk_wind.gs_indx , fontsize=10, color = 'k')

                        # plt.text(dlnk_wind_start+0.15, bottom_vert_loc+0.1, dlnk_wind.window_ID , fontsize=10, color = 'k')

                        # save off the rotator choice so that we can look it up again
                        dlnk_choices_rectangle_rotator_hist[dlnk_wind] = dlnk_choices_rectangle_rotator

                        dlnk_choices_rectangle_rotator =  (dlnk_choices_rectangle_rotator+2)%dlnk_rotation_rollover

                        num_dlnk_wind += 1

            # plot the executed down links
            if self.winds_plot_dlnks:
                num_dlnk_exe = 0
                if len(dlnk_winds_flat) > 0:
                    for dlnk_wind in dlnk_winds_flat[sat_indx]:

                        dlnk_start = (dlnk_wind.start-plot_start).total_seconds()/time_divisor
                        dlnk_end = (dlnk_wind.end-plot_start).total_seconds()/time_divisor
                        middle = (dlnk_end + dlnk_start)/2

                        gs_indx = dlnk_wind.gs_indx

                        #  update the rotator value if we've already added this window to the plot in the "choices" code above
                        if dlnk_wind in dlnk_choices_rectangle_rotator_hist.keys ():
                            dlnk_rectangle_rotator = dlnk_choices_rectangle_rotator_hist[dlnk_wind]

                        # plot the task duration
                        bottom_vert_loc = dlnk_rectangle_rotator
                        radius =  dv_radius_scale(obs_wind.scheduled_data_vol)
                        d = Circle((middle, bottom_vert_loc), radius,alpha=1,fill=False,color='#0000FF',hatch='///////')
                        current_axis.add_patch(d)

                        dlnk_rectangle_rotator =  (dlnk_rectangle_rotator+2)%dlnk_rotation_rollover

                        if dlnk_wind.window_ID in all_wind_ids:
                            raise Exception('Found a duplicate unique window ID where it should not have been possible')
                        all_wind_ids.append(dlnk_wind.window_ID)

                        if plot_include_labels:
                            label_text = ""

                            # if we have been given route indices, then include them in the label
                            if (len(route_indcs_by_wind.keys())) > 0:
                                dr_indcs = route_indcs_by_wind[dlnk_wind]
                                for dr_indx in dr_indcs:
                                    label_text += "%d,"%(dr_indx)
                                label_text += ";"
                            
                            #  add the ground station index to the label
                            label_text += "%d"%(gs_indx)

                            #   put label in desired vertical spot
                            left_horizontal_loc = dlnk_start + 0.15
                            
                            if dlnk_label_rotator == 0:
                                plt.text(left_horizontal_loc, 0.1, label_text , fontsize=10, color = 'k')
                            elif dlnk_label_rotator == 1:
                                plt.text( left_horizontal_loc, 1.1, label_text , fontsize=10, color = 'k')

                            dlnk_label_rotator = (dlnk_label_rotator+1)%dlnk_rotation_rollover

                        num_dlnk_exe += 1

            ###################
            # xlnks
            ###################

            #  plot the potential cross-links
            if self.winds_plot_xlnks_choices:
                num_xlnk_wind = 0

                if len(all_xlnk_winds_flat) > 0:
                    for xlnk_wind in all_xlnk_winds_flat[sat_indx]:

                        xlnk_wind_start = (xlnk_wind.start-plot_start).total_seconds()/time_divisor
                        xlnk_wind_end = (xlnk_wind.end-plot_start).total_seconds()/time_divisor
                        middle = (xlnk_wind_end + xlnk_wind_start)/2

                        # plot the task duration
                        bottom_vert_loc = xlnk_choices_rectangle_rotator
                        radius =  dv_radius_scale(xlnk_wind.data_vol)
                        x_w = Circle((middle, bottom_vert_loc), radius,alpha=1,fill=True,color='#FFBCBC')

                        current_axis.add_patch(x_w)

                        # plt.text(xlnk_wind_start+0.15, bottom_vert_loc+0.1, xlnk_wind.window_ID , fontsize=10, color = 'k')

                        # save off the rotator choice so that we can look it up again
                        xlnk_choices_rectangle_rotator_hist[xlnk_wind] = xlnk_choices_rectangle_rotator

                        xlnk_choices_rectangle_rotator =  (xlnk_choices_rectangle_rotator+2)%xlnk_rotation_rollover
                        num_xlnk_wind += 1

            #  plot the executed cross-links
            if self.winds_plot_xlnks:
                num_xlnk_exe = 0
                if len(xlnk_winds_flat) > 0:
                    for xlnk_wind in xlnk_winds_flat[sat_indx]:

                        xlnk_start = (xlnk_wind.start-plot_start).total_seconds()/time_divisor
                        xlnk_end = (xlnk_wind.end-plot_start).total_seconds()/time_divisor
                        middle = (xlnk_end + xlnk_start)/2

                        #  update the rotator value if we've already added this window to the plot in the "choices" code above
                        if xlnk_wind in xlnk_choices_rectangle_rotator_hist.keys ():
                            xlnk_rectangle_rotator = xlnk_choices_rectangle_rotator_hist[xlnk_wind]

                        bottom_vert_loc = xlnk_rectangle_rotator

                        # if we have been given route indices, then then use one of them to configure the crosslink box color
                        dr_indx = None
                        if (len(route_indcs_by_wind.keys())) > 0:
                            dr_indx = route_indcs_by_wind[xlnk_wind][self.route_index_to_use]
                            xlnk_color_indx = dr_indx %  self.xlnk_color_rollover
                        #  otherwise just go with the first color
                        else:
                            xlnk_color_indx = 0

                        xlnk_color = self.xlnk_colors[xlnk_color_indx]
                        
                        # plot the task duration
                        radius =  dv_radius_scale(xlnk_wind.scheduled_data_vol)
                        x = Circle((middle, bottom_vert_loc), radius,alpha=1,fill=False,color=xlnk_color,hatch='///////')
                        current_axis.add_patch(x)

                        xlnk_rectangle_rotator =  (xlnk_rectangle_rotator+2)%xlnk_rotation_rollover

                        if plot_include_labels:
                            other_sat_indx = xlnk_wind.xsat_indx if xlnk_wind.xsat_indx != sat_indx else xlnk_wind.sat_indx

                            #   put label in desired vertical spot
                            left_horizontal_loc = xlnk_start + 0.15

                            #  again, if we know route indices include them in label
                            if dr_indx:
                                label_text = "%d,%d" %(dr_indx,other_sat_indx)
                            else:                                
                                label_text = "%d" %(other_sat_indx)

                            if xlnk_label_rotator == 0:
                                plt.text(left_horizontal_loc, 0.1, label_text , fontsize=10, color = 'k')
                            elif xlnk_label_rotator == 1:
                                plt.text( left_horizontal_loc, 1.1, label_text , fontsize=10, color = 'k')

                            xlnk_label_rotator = (xlnk_label_rotator+1)%xlnk_rotation_rollover

            #  if were at the last satellite ( at the bottom of all the plots), then add X axis labels
            if not plot_indx+1 == num_sats:
                ax = plt.gca()
                plt.setp(ax.get_xticklabels(), visible=False)

        legend_objects = []
        legend_objects_labels = []
        if d_w: 
            legend_objects.append(d_w)
            legend_objects_labels.append('D all')
        if d: 
            legend_objects.append(d)
            legend_objects_labels.append('Dlnk')
        if x_w: 
            legend_objects.append(x_w)
            legend_objects_labels.append('X all')
        if x: 
            legend_objects.append(x)
            legend_objects_labels.append('Xlnk')

        plt.legend(legend_objects, legend_objects_labels ,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.xlabel('Time (%s)'%(self.time_units))

        # use the last axes to set the entire plot background color
        axes.patch.set_facecolor('w')

        if show:
            plt.show()
        else:
            savefig(fig_name,format=self.plot_fig_extension)

    def plot_energy_usage(
        self,
        sats_ids_list,
        energy_usage,
        ecl_winds,
        plot_start,
        plot_end,
        plot_title = 'Energy Utilization', 
        plot_size_inches = (12,12),
        show=False,
        fig_name='plots/energy_plot.pdf'):

        plot_labels = {
            "e usage": "e usage",
            "e max": "e max",
            "e min": "e min",
            "ecl": "ecl"
        }

        if self.time_units == 'hours':
            time_divisor = 3600
        if self.time_units == 'minutes':
            time_divisor = 60
        
        time_to_end = (plot_end-plot_start).total_seconds()/time_divisor

        num_sats = len(sats_ids_list)

        #  make a new figure
        plt.figure()

        #  create subplots for satellites
        fig = plt.gcf()
        fig.set_size_inches( plot_size_inches)
        # print fig.get_size_inches()

        plt.title( plot_title)

        #  these hold the very last plot object of a given type added. Used for legend below
        e_usage_plot = None
        e_max_plot = None
        e_min_plot = None
        ecl_plot = None

        # for each agent
        for  plot_indx, sat_id in enumerate (sats_ids_list):
            #  get the index for this ID
            sat_indx = self.sat_id_order.index(str(sat_id))

            #  make a subplot for each
            axes = plt.subplot( num_sats,1,plot_indx+1)
            if plot_indx == np.floor(num_sats/2):
                plt.ylabel('Satellite Index\n\n' + str(sat_indx))
            else:
                plt.ylabel('' + str(sat_indx))


            # no y-axis labels
            plt.tick_params(
                axis='y',
                which='both',
                left='off',
                right='off',
                labelleft='off'
            )

            # set axis length. Time starts at 0
            vert_min = self.energy_usage_plot_params['plot_bound_e_min_Wh']
            vert_max = self.energy_usage_plot_params['plot_bound_e_max_Wh']
            plt.axis((0, time_to_end, vert_min, vert_max))

            current_axis = plt.gca()

            # the first return value is a handle for our line, everything else can be ignored
            if energy_usage:
                e_usage_plot,*dummy = plt.plot(energy_usage['time_mins'],energy_usage['e_sats'][sat_indx], label =  plot_labels["e usage"])

            if self.energy_usage_plot_params['include_eclipse_windows']:
                for ecl_wind in ecl_winds[sat_indx]:
                    ecl_wind_start = (ecl_wind.start-plot_start).total_seconds()/time_divisor
                    ecl_wind_end = (ecl_wind.end-plot_start).total_seconds()/time_divisor

                    height = 2
                    ecl_plot = Rectangle((ecl_wind_start, vert_min), ecl_wind_end-ecl_wind_start, vert_min+height,alpha=0.3,fill=True,color='#202020',label= plot_labels["ecl"])

                    current_axis.add_patch(ecl_plot)


        legend_objects = []
        if e_usage_plot: 
            legend_objects.append(e_usage_plot)
        if e_max_plot: 
            legend_objects.append(e_max_plot)
        if e_min_plot: 
            legend_objects.append(e_min_plot)
        if ecl_plot: 
            legend_objects.append(ecl_plot)

        plt.legend(handles=legend_objects ,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.xlabel('Time (%s)'%(self.time_units))

        # use the last axes to set the entire plot background color
        axes.patch.set_facecolor('w')

        if show:
            plt.show()
        else:
            savefig(fig_name,format=self.plot_fig_extension)


    def plot_data_usage(
        self,
        sats_ids_list,
        data_usage,
        ecl_winds,
        plot_start,
        plot_end,
        plot_title = 'data Utilization', 
        plot_size_inches = (12,12),
        show=False,
        fig_name='plots/data_plot.pdf'):

        plot_labels = {
            "d usage": "d usage",
            "d max": "d max",
            "d min": "d min",
            "ecl": "ecl"
        }

        if self.time_units == 'hours':
            time_divisor = 3600
        if self.time_units == 'minutes':
            time_divisor = 60
        
        time_to_end = (plot_end-plot_start).total_seconds()/time_divisor

        num_sats = len(sats_ids_list)

        #  make a new figure
        plt.figure()

        #  create subplots for satellites
        fig = plt.gcf()
        fig.set_size_inches( plot_size_inches)
        # print fig.get_size_inches()

        plt.title( plot_title)

        #  these hold the very last plot object of a given type added. Used for legend below
        d_usage_plot = None
        d_max_plot = None
        d_min_plot = None
        ecl_plot = None

        # for each agent
        for  plot_indx, sat_id in enumerate (sats_ids_list):
            #  get the index for this ID
            sat_indx = self.sat_id_order.index(str(sat_id))

            #  make a subplot for each
            axes = plt.subplot( num_sats,1,plot_indx+1)
            if plot_indx == np.floor(num_sats/2):
                plt.ylabel('Satellite Index\n\n' + str(sat_indx))
            else:
                plt.ylabel('' + str(sat_indx))


            # no y-axis labels
            plt.tick_params(
                axis='y',
                which='both',
                left='off',
                right='off',
                labelleft='off'
            )

            # set axis length. Time starts at 0
            vert_min = self.data_usage_plot_params['plot_bound_d_min_Gbit']
            vert_max = self.data_usage_plot_params['plot_bound_d_max_Gbit']
            plt.axis((0, time_to_end, vert_min, vert_max))

            current_axis = plt.gca()

            # the first return value is a handle for our line, everything else can be ignored
            if data_usage:
                # data_usage['d_sats'][sat_indx] is currently in Mb...
                sat_data_usage_Gb = [num/1000 for num in data_usage['d_sats'][sat_indx]]
                d_usage_plot,*dummy = plt.plot(data_usage['time_mins'],sat_data_usage_Gb, label =  plot_labels["d usage"])

            if self.data_usage_plot_params['include_eclipse_windows']:
                for ecl_wind in ecl_winds[sat_indx]:
                    ecl_wind_start = (ecl_wind.start-plot_start).total_seconds()/time_divisor
                    ecl_wind_end = (ecl_wind.end-plot_start).total_seconds()/time_divisor

                    height = 2
                    ecl_plot = Rectangle((ecl_wind_start, vert_min), ecl_wind_end-ecl_wind_start, vert_min+height,alpha=0.3,fill=True,color='#202020',label= plot_labels["ecl"])

                    current_axis.add_patch(ecl_plot)


        legend_objects = []
        if d_usage_plot: 
            legend_objects.append(d_usage_plot)
        if d_max_plot: 
            legend_objects.append(d_max_plot)
        if d_min_plot: 
            legend_objects.append(d_min_plot)
        if ecl_plot: 
            legend_objects.append(ecl_plot)

        plt.legend(handles=legend_objects ,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.xlabel('Time (%s)'%(self.time_units))

        # use the last axes to set the entire plot background color
        axes.patch.set_facecolor('w')

        if show:
            plt.show()
        else:
            savefig(fig_name,format=self.plot_fig_extension)

    def plot_obs_aoi(
        self,
        targ_ids_list,
        aoi_curves_by_targID,
        plot_start,
        plot_end,
        plot_title = 'Observation Target AoI', 
        plot_size_inches = (12,12),
        show=False,
        fig_name='plots/obs_aoi_plot.pdf'):

        plot_labels = {
            "aoi": "aoi",
        }

        if self.time_units == 'hours':
            time_divisor = 3600
        if self.time_units == 'minutes':
            time_divisor = 60
        
        time_to_end = (plot_end-plot_start).total_seconds()/time_divisor

        num_targs = len(targ_ids_list)

        #  make a new figure
        plt.figure()

        #  create subplots for satellites
        fig = plt.gcf()
        fig.set_size_inches( plot_size_inches)
        # print fig.get_size_inches()

        plt.title( plot_title)

        #  these hold the very last plot object of a given type added. Used for legend below
        aoi_plot = None

        # for each agent
        for  plot_indx, targ_id in enumerate (targ_ids_list):
            #  get the index for this ID
            targ_indx = self.all_targ_IDs.index(targ_id)

            #  make a subplot for each
            axes = plt.subplot( num_targs,1,plot_indx+1)
            if plot_indx == np.floor(num_targs/2):
                plt.ylabel('Target Index\n\n' + str(targ_indx))
            else:
                plt.ylabel('' + str(targ_indx))


            # no y-axis labels
            plt.tick_params(
                axis='y',
                which='both',
                left='off',
                right='off',
                labelleft='off'
            )

            # set axis length. Time starts at 0
            vert_min = self.obs_aoi_metrics_plot_params['plot_bound_min_aoi_h']
            vert_max = self.obs_aoi_metrics_plot_params['plot_bound_max_aoi_h']
            plt.axis((0, time_to_end, vert_min, vert_max))

            current_axis = plt.gca()

            # the first return value is a handle for our line, everything else can be ignored
            aoi_plot,*dummy = plt.plot(aoi_curves_by_targID[targ_id]['x'],aoi_curves_by_targID[targ_id]['y'], label =  plot_labels["aoi"])

        legend_objects = []
        if aoi_plot: 
            legend_objects.append(aoi_plot)

        plt.legend(handles=legend_objects ,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.xlabel('Time (%s)'%(self.time_units))

        # use the last axes to set the entire plot background color
        axes.patch.set_facecolor('w')

        if show:
            plt.show()
        else:
            savefig(fig_name,format=self.plot_fig_extension)

    def plot_sat_tlm_cmd_aoi(
        self,
        sats_ids_list,
        aoi_curves_by_sat_indx,
        aoi_option,
        plot_start,
        plot_end,
        plot_title = 'Satellite TLM/CMD Down/Uplink AoI', 
        plot_size_inches = (12,12),
        show=False,
        fig_name='plots/sat_aoi_plot.pdf'):

        plot_labels = {
            "aoi": "aoi",
        }

        if self.time_units == 'hours':
            time_divisor = 3600
        if self.time_units == 'minutes':
            time_divisor = 60
        
        time_to_end = (plot_end-plot_start).total_seconds()/time_divisor

        num_sats = len(sats_ids_list)

        #  make a new figure
        plt.figure()

        #  create subplots for satellites
        fig = plt.gcf()
        fig.set_size_inches( plot_size_inches)
        # print fig.get_size_inches()

        plt.title( plot_title)

        #  these hold the very last plot object of a given type added. Used for legend below
        aoi_plot = None

        # for each agent
        for  plot_indx, sat_id in enumerate (sats_ids_list):
            sat_indx = self.sat_id_order.index(str(sat_id))


            #  make a subplot for each
            axes = plt.subplot( num_sats,1,plot_indx+1)
            if plot_indx == np.floor(num_sats/2):
                plt.ylabel('Satellite Index\n\n' + str(sat_indx))
            else:
                plt.ylabel('' + str(sat_indx))


            # no y-axis labels
            plt.tick_params(
                axis='y',
                which='both',
                left='off',
                right='off',
                labelleft='off'
            )

            # set axis length. Time starts at 0
            if aoi_option == 'cmd':
                vert_min = self.cmd_aoi_metrics_plot_params['plot_bound_min_aoi_h']
                vert_max = self.cmd_aoi_metrics_plot_params['plot_bound_max_aoi_h']
            elif aoi_option == 'tlm':
                vert_min = self.tlm_aoi_metrics_plot_params['plot_bound_min_aoi_h']
                vert_max = self.tlm_aoi_metrics_plot_params['plot_bound_max_aoi_h']
            plt.axis((0, time_to_end, vert_min, vert_max))

            current_axis = plt.gca()

            # the first return value is a handle for our line, everything else can be ignored
            aoi_plot,*dummy = plt.plot(aoi_curves_by_sat_indx[sat_indx]['x'],aoi_curves_by_sat_indx[sat_indx]['y'], label =  plot_labels["aoi"])

        legend_objects = []
        if aoi_plot: 
            legend_objects.append(aoi_plot)

        plt.legend(handles=legend_objects ,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.xlabel('Time (%s)'%(self.time_units))

        # use the last axes to set the entire plot background color
        axes.patch.set_facecolor('w')

        if show:
            plt.show()
        else:
            savefig(fig_name,format=self.plot_fig_extension)

    @staticmethod
    def plot_route_latdv_pareto(all_routes,weights_tups,plot_name,show= False):

        # print ('obs.data_vol, total dlnk dv, ratio')
        # print (obs.data_vol,sum(dr.data_vol for dr in routes),sum(dr.data_vol for dr in routes)/obs.data_vol)
        # print ('min latency, ave latency, max latency')
        all_ave_latencies = []
        all_cum_dv = []
        for routes in all_routes:
            raise RuntimeWarning("TODO:  update the line of code below")
            latencies = [(rt.route[-1].end - rt.route[0].end).total_seconds()/60 for rt in  routes]
            all_ave_latencies.append ( np.mean(latencies))
            all_cum_dv.append ( sum(dr.data_vol for dr in routes))
        # print (np.min( latencies),np.mean( latencies),np.max( latencies))

        #  make a new figure
        plt.figure()
        plt.plot(all_ave_latencies,all_cum_dv,'o')

        labels = ['{0}'.format(i) for i in range(len(all_ave_latencies))]
        # labels = ['{0} dv{1} la{2}'.format(i,tup[0],tup[2]) for i, tup in   enumerate(weights_tups)]
        indx = 0
        for label, x, y in zip(labels, all_ave_latencies, all_cum_dv):
            
            if indx%2 ==0: pos = (30, 0)
            else: pos = (30, -20)

            if indx == 1: pos=(30,0)
            if indx == 2: pos=(30,0)
            if indx == 3: pos=(-20,0)
            if indx == 4: pos=(-20,-20)
            if indx == 5: pos=(-20,-40)
            if indx == 6: pos=(-20,-60)
            if indx == 7: pos=(30,-20)
            if indx == 8: pos=(-20,-80)
            if indx == 9: pos=(30,-40)

            plt.annotate(
                label,
                xy=(x, y), xytext=pos,
                textcoords='offset points', ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
            indx += 1

        fig = plt.gcf()
        fig.set_size_inches(12,8)
        plt.xlabel('Average latency (mins)')
        plt.ylabel('Total DV (Mb)')
        plt.title( 'Ave Latency - Total DV (over all paths) Pareto frontier')
        if show:
            plt.show()
        else:
            savefig(plot_name,format='pdf')