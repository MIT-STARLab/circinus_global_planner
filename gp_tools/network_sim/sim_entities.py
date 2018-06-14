# Objects used in network simulation
#
# @author: Kit Kennedy

from collections import namedtuple
import pdb

import circinus_tools.metrics.metrics_utils as met_util

class NetSimEntity():
    def __init__(self,num_sats,num_gs,start_dt,end_dt):

        self.time_units = 'seconds'
        self.start_dt = start_dt
        self.end_dt = end_dt
        # start at zero seconds ( arbitrary, but works well)
        self.start_time_s = 0
        start = self.start_time_s
        self.end_time_s = (end_dt-start_dt).total_seconds ()
        self.sat_update_hist = [met_util.UpdateHistory(t=[start],last_update_time=[start]) for i in range(num_sats)]
        self.gs_update_hist = [met_util.UpdateHistory(t=[start],last_update_time=[start]) for i in range(num_gs)]

    def pull_update_times(self, current_time,other,option='sat_times'):
        if option =='sat_times':
            update_hist = self.sat_update_hist
            other_update_hist =  other.sat_update_hist
        elif option =='gs_times':
            update_hist = self.gs_update_hist
            other_update_hist =  other.gs_update_hist
        else:
            raise NotImplementedError

        for indx,uh in enumerate (update_hist):
            # print ('yo')
            # print (current_time)
            # print (update_hist[indx])
            # print (self.sat_update_hist)
            # pdb.set_trace()
            last_ut =  update_hist[indx].last_update_time[-1]
            other_last_ut =  other_update_hist[indx].last_update_time[-1]
            
            #  if our last update time is before the others last update time for this index, then we should add that more recent update time to our history
            #  note that this should never end up updating the history for self's sat or gs indx
            if  last_ut <  other_last_ut:
                update_hist[indx].t.append ( current_time)
                update_hist[indx].last_update_time.append (other_last_ut)
                # print ('update_hist[indx].t')
                # print (update_hist[indx].t)

    def pull_ttc(self, current_time, other):
        """ pull command and telemetry data update times from other entity
        
        updates all of the update times on the entity doing the pulling (e.g. during a cross-link or downlink).  note that the update_self method should be called on all relevant parties in advance ( which may be more than just the two entities involved in this specific exchange)

        note that this may be used to pull TTC from a satellite to a ground station, in which case it might not make that much sense to update all of the gs_times. but whatever, shouldn't hurt.

        :param other:  the other entity ( satellite or ground station)
        :type other: [NetSimEntity]
        """
        self.pull_update_times(current_time,other,option ='sat_times')
        self.pull_update_times(current_time,other,option ='gs_times')

    def  refresh_update_time( self,units='seconds',time_option ='relative_to_start'):

        if time_option == 'relative_to_start':
            pass
        # TODO: should implement this option - helpful for case where you want to get datetimes back
        elif time_option == 'absolute':
            raise NotImplementedError
        else:
            raise NotImplementedError

class NetSimSat(NetSimEntity):
    def __init__(self,sat_indx,sat_id,num_sats,num_gs,start_dt,end_dt):
        self.sat_indx = sat_indx
        self.ID = sat_id
        super(NetSimSat, self).__init__(num_sats,num_gs,start_dt,end_dt)

    def  refresh_update_time( self,current_time,units='seconds',time_option ='relative_to_start'):
        #  don't append to the update time list for self -  just replace the first element.  otherwise we'd have a list as long as the number of time elements
        self.sat_update_hist[self.sat_indx].last_update_time[0] = current_time
        super(). refresh_update_time(units,time_option)

    def get_merged_cmd_update_hist(self,gs_agents,gs_id_ignore_list):
        """ gets the merged command update history for this satellite
        
        Gets the update history for every ground station as seen by this satellite, and then merges these into a single update history. The update history for ground station gs_indx for this satellite  is a recording of when this satellite last heard from that ground station. By merging across all ground stations, we get a recording of when the satellite last heard from any ground station, which we assume is a good proxy for when ground commanding was last updated. ( note the underlying assumption that all ground stations have equal relevance for commanding the satellite)
        :param all_gs_IDs: [description]
        :type all_gs_IDs: [type]
        :param gs_id_ignore_list: [description]
        :type gs_id_ignore_list: [type]
        :returns: [description]
        :rtype: {[type]}
        """

        update_hists = []

        # grab the update time histories from all the ground stations
        for gs_agent in gs_agents:
            # if we're ignoring this ground station, then continue
            if gs_agent.ID in gs_id_ignore_list:
                continue

            update_hists.append (self.gs_update_hist[gs_agent.gs_indx])

        return met_util.merge_update_histories ( update_hists,self.end_time_s)

class NetSimGS(NetSimEntity):
    def __init__(self,gs_indx,gs_id,num_sats,num_gs,start_dt,end_dt):
        self.gs_indx = gs_indx
        self.ID = gs_id
        super(NetSimGS, self).__init__(num_sats,num_gs,start_dt,end_dt)

    def  refresh_update_time( self,current_time,units='seconds',time_option ='relative_to_start'):
        #  don't append to the update time list for self -  just replace the first element.  otherwise we'd have a list as long as the number of time elements
        self.gs_update_hist[self.gs_indx].last_update_time[0] = current_time
        super(). refresh_update_time(units,time_option)


    def get_sat_tlm_update_hist(self,sat_agent):
        """gets the telemetry update history for a single satellite
        
        Gets the update history for a single satellite as seen by this ground station. The update history for satellite sat_indx for this ground station is a recording of when this ground station last heard from that satellite. By merging the update history returned here across satellites, you can get the merged telemetry update history for a single satellite across the full ground station network
        :param sat_indx: [description]
        :type sat_indx: [type]
        :returns: [description]
        :rtype: {[type]}
        """

        return self.sat_update_hist[sat_agent.sat_indx]
        
        



