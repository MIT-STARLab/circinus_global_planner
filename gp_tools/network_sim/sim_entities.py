# Objects used in network simulation
#
# @author: Kit Kennedy

from collections import namedtuple
import pdb

UpdateHistory = namedtuple('UpdateHistory', 't last_update_time')

class NetSimEntity():
    def __init__(self,num_sats,num_gs,start_dt,end_dt):

        self.time_units = 'seconds'
        self.start_dt = start_dt
        self.end_dt = end_dt
        # start at zero seconds ( arbitrary, but works well)
        self.start_time_s = 0
        start = self.start_time_s
        self.end_time_s = (end_dt-start_dt).total_seconds ()
        self.sat_update_hist = [UpdateHistory(t=[start],last_update_time=[start]) for i in range(num_sats)]
        self.gs_update_hist = [UpdateHistory(t=[start],last_update_time=[start]) for i in range(num_gs)]

    @staticmethod
    def fix_update_times(t_lut_zip_sorted,last_time):
        """[summary]
        
        this algorithm is the same as that in gp_metrics.py, get_av_aoi_routing(). "last_update_time" is the most recent time that the originating entity updated itself (i.e.  it's a creation time) and "t" is the time at which this last update time was delivered to self (i.e. it's a delivery time). We're producing a delivery/creation matrix with this code

        todo: this is suboptimal code reuse.  should merge with the stuff in GP metrics
        :param t_lut_zip_sorted: [description]
        :type t_lut_zip_sorted: [type]
        :param last_time: [description]
        :type last_time: [type]
        """

        #  start things off with the very first point
        merged_t = [t_lut_zip_sorted[0][0]]
        merged_lut = [t_lut_zip_sorted[0][1]]

        curr_time = merged_t[-1]
        curr_lut = merged_lut[-1]
        for next_item in t_lut_zip_sorted:
            next_time = next_item[0]
            next_lut = next_item[1]
            # if we received a more recent last update time at the current time, then we need to fix the last update time in the merged list
            if next_time == curr_time:
                if next_lut > curr_lut:
                    curr_lut = next_lut
                    merged_lut[-1] = curr_lut

            #  if the next point is a later time and also has a more recent last update time, then add to merged
            if next_time > curr_time:
                if next_lut > curr_lut:
                    curr_time = next_time
                    curr_lut = next_lut
                    merged_t.append(curr_time)
                    merged_lut.append(curr_lut)

        # want to bookend this with a last time so it's explicit what time window we're looking at
        if last_time > curr_time:
            merged_t.append(last_time)
            merged_lut.append(curr_lut)            

        return merged_t,merged_lut

    @staticmethod
    def  merge_update_histories ( update_hists, end_time_s):
        merged_t_dirty = []
        # lut = last update time
        merged_lut_dirty = []

        #  first merge the update time histories
        for update_hist in update_hists:
            merged_t_dirty += update_hist.t
            merged_lut_dirty += update_hist.last_update_time

        # sort by t value
        zip_sorted = sorted(zip(merged_t_dirty,merged_lut_dirty),key= lambda z: z[0])

        #  there could still be duplicate t values, so let's get rid of those. Want to grab the best (most recent) lut value at each t
        merged_t,merged_lut =  NetSimEntity.fix_update_times(zip_sorted, end_time_s)
        
        return UpdateHistory(t=merged_t,last_update_time=merged_lut)

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
    def __init__(self,sat_indx,num_sats,num_gs,start_dt,end_dt):
        self.sat_indx = sat_indx
        super(NetSimSat, self).__init__(num_sats,num_gs,start_dt,end_dt)

    def  refresh_update_time( self,current_time,units='seconds',time_option ='relative_to_start'):
        #  don't append to the update time list for self -  just replace the first element.  otherwise we'd have a list as long as the number of time elements
        self.sat_update_hist[self.sat_indx].last_update_time[0] = current_time
        super(). refresh_update_time(units,time_option)

    def get_merged_cmd_update_hist(self,num_gs,all_gs_IDs,gs_id_ignore_list):
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
        for gs_indx in  range (num_gs):
            # if we're ignoring this ground station, then continue
            if all_gs_IDs[gs_indx] in gs_id_ignore_list:
                continue

            update_hists.append (self.gs_update_hist[gs_indx])

        return self.merge_update_histories ( update_hists,self.end_time_s)

class NetSimGS(NetSimEntity):
    def __init__(self,gs_indx,num_sats,num_gs,start_dt,end_dt):
        self.gs_indx = gs_indx
        super(NetSimGS, self).__init__(num_sats,num_gs,start_dt,end_dt)

    def  refresh_update_time( self,current_time,units='seconds',time_option ='relative_to_start'):
        #  don't append to the update time list for self -  just replace the first element.  otherwise we'd have a list as long as the number of time elements
        self.gs_update_hist[self.gs_indx].last_update_time[0] = current_time
        super(). refresh_update_time(units,time_option)


    def get_sat_tlm_update_hist(self,sat_indx):
        """gets the telemetry update history for a single satellite
        
        Gets the update history for a single satellite as seen by this ground station. The update history for satellite sat_indx for this ground station is a recording of when this ground station last heard from that satellite. By merging the update history returned here across satellites, you can get the merged telemetry update history for a single satellite across the full ground station network
        :param sat_indx: [description]
        :type sat_indx: [type]
        :returns: [description]
        :rtype: {[type]}
        """

        return self.sat_update_hist[sat_indx]
        
        



