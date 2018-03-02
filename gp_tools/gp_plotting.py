#  plotting tools for global planner
# 
# @author Kit Kennedy

import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig
from matplotlib.patches import Rectangle
import numpy as np

from .routing_objects import DataRoute

class GPPlotting():
    """docstring for GP route selection"""
    def __init__(self,params):

        self.plot_fig_extension=params['plot_fig_extension']
        self.plot_include_labels=params['plot_include_labels']
        self.time_units=params['time_units']
        self.links_plot_dlnks=params['links_plot_dlnks']
        self.links_plot_dlnks_choices=params['links_plot_dlnks_choices']
        self.links_plot_xlnks=params['links_plot_xlnks']
        self.links_plot_xlnks_choices=params['links_plot_xlnks_choices']
        
        # self.num_sats=params['num_sats']
        # self.num_paths=params['route_selection_num_paths']
        # self.start_utc_dt  =params['start_utc_dt']
        # self.end_utc_dt  =params['end_utc_dt']

    def plot_links(
        self,
        num_sats,
        all_dlnk_winds_flat,
        dlnk_winds_flat, 
        all_xlnk_winds_flat,
        xlnk_winds_flat,
        plot_start,
        plot_end,
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

        #  create subplots for satellites
        axes = plt.subplot(num_sats,1,1)
        axes.patch.set_facecolor('w')
        fig = plt.gcf()
        fig.set_size_inches((12,6))
        # print fig.get_size_inches()

        plt.title('Links Plot')

        # have to sort before applying axis labels, otherwise x label shows up in a weird place
        # sats.sort(key=lambda x: x.agent_ID)

        # for each agent
        for sat_indx in range (num_sats):

            # 
            plt.subplot( num_sats,1,sat_indx+1)
            if sat_indx == np.floor(num_sats/2):
                plt.ylabel('Satellite Number\n\n' + str(sat_indx) + "                ")
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

            #  this is used to alternate the vertical position of labels
            label_rotator = 0

            #  these hold the very last plot object of a given type added. Used for legend below
            d_w = None
            d = None
            x_w = None
            x = None

            ###################
            # dlnks
            ###################

            #  plot the potential down links
            if self.links_plot_dlnks_choices:
                num_dlnk_wind = 0

                if len(all_dlnk_winds_flat) > 0:
                    for dlnk_wind in all_dlnk_winds_flat[sat_indx]:

                        dlnk_wind_start = (dlnk_wind.start-plot_start).total_seconds()/time_divisor
                        dlnk_wind_end = (dlnk_wind.end-plot_start).total_seconds()/time_divisor

                        # plot the task duration
                        d_w = Rectangle((dlnk_wind_start, 0), dlnk_wind_end-dlnk_wind_start, 2,alpha=1,fill=True,color='#BFBFFF')

                        current_axis.add_patch(d_w)
                        plt.text( (dlnk_wind_end+dlnk_wind_start)/2 - 0.15, 0.1, dlnk_wind.gs_ID , fontsize=10, color = 'k')

                        num_dlnk_wind += 1

            # plot the executed down links
            if self.links_plot_dlnks:
                num_dlnk_exe = 0
                if len(dlnk_winds_flat) > 0:
                    for dlnk_wind in dlnk_winds_flat[sat_indx]:

                        dlnk_start = (dlnk_wind.start-plot_start).total_seconds()/time_divisor
                        dlnk_end = (dlnk_wind.end-plot_start).total_seconds()/time_divisor

                        # plot the task duration
                        d = Rectangle((dlnk_start, 0), dlnk_end-dlnk_start, 2,alpha=1,fill=False,color='#0000FF',hatch='///////')
                        current_axis.add_patch(d)

                        num_dlnk_exe += 1

            ###################
            # xlnks
            ###################

            #  plot the potential cross-links
            if self.links_plot_xlnks_choices:
                num_xlnk_wind = 0

                if len(all_xlnk_winds_flat) > 0:
                    for xlnk_wind in all_xlnk_winds_flat[sat_indx]:

                        xlnk_wind_start = (xlnk_wind.start-plot_start).total_seconds()/time_divisor
                        xlnk_wind_end = (xlnk_wind.end-plot_start).total_seconds()/time_divisor

                        # plot the task duration
                        x_w = Rectangle((xlnk_wind_start, 0), xlnk_wind_end-xlnk_wind_start, 2,alpha=1,fill=True,color='#FFBCBC')

                        current_axis.add_patch(x_w)

                        num_xlnk_wind += 1

            #  plot the executed cross-links
            if self.links_plot_xlnks:
                num_xlnk_exe = 0
                if len(xlnk_winds_flat) > 0:
                    for xlnk_wind in xlnk_winds_flat[sat_indx]:

                        xlnk_start = (xlnk_wind.start-plot_start).total_seconds()/time_divisor
                        xlnk_end = (xlnk_wind.end-plot_start).total_seconds()/time_divisor

                        # plot the task duration
                        x = Rectangle((xlnk_start, 0), xlnk_end-xlnk_start, 2,alpha=1,fill=False,color='#FF0000',hatch='///////')
                        current_axis.add_patch(x)

                        if self.plot_include_labels:
                            other_sat_indx = xlnk_wind.xsat_indx if xlnk_wind.xsat_indx != sat_indx else xlnk_wind.sat_indx

                            #   put label in desired vertical spot
                            if label_rotator == 0:
                                plt.text( (xlnk_end+xlnk_start)/2 - 0.15, 0.1, str(other_sat_indx) , fontsize=10, color = 'k')
                            elif label_rotator == 1:
                                plt.text( (xlnk_end+xlnk_start)/2 - 0.15, 0.6, str(other_sat_indx) , fontsize=10, color = 'k')
                            elif label_rotator == 2:
                                plt.text( (xlnk_end+xlnk_start)/2 - 0.15, 1.1, str(other_sat_indx) , fontsize=10, color = 'k')
                            else:
                                plt.text( (xlnk_end+xlnk_start)/2 - 0.15, 1.6, str(other_sat_indx) , fontsize=10, color = 'k')

                            label_rotator = (label_rotator+1)%4

            #  if were at the last satellite ( at the bottom of all the plots), then add X axis labels
            if not sat_indx+1 == num_sats:
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

        if show:
            plt.show()
        else:
            savefig(fig_name,format=self.plot_fig_extension)