__author__="huziy"
__date__ ="$Sep 17, 2011 5:12:55 PM$"

import os
import members
import data_select
import gevfit
from datetime import datetime
from datetime import timedelta


import pickle
import numpy as np

def gev_fit_all_members(high_flow = True, member_ids = [], 
                        data_folder = '',
                        file_name_pattern = '',
                        start_date = None,
                        end_date = None,
                        start_month = 1,
                        end_month = 12,
                        duration_days = timedelta(days = 1)):
    """
    gev fit using data from all members
    data_folder - path to the folder with input data (streamflow)
    start_month -
    end_month - 
    """

    param_file = 'high' if high_flow else 'low'
    for id in member_ids:
        param_file += '_' + id
    if os.path.isfile(param_file):
        print('delete {0}, to reoptimize'.format(param_file))
        return pickle.load(open(param_file))

    #select data
    path_pattern = os.path.join(data_folder, file_name_pattern)
    all_extremes = []
    for id in member_ids:
        print(id)
        the_path = path_pattern.format(id)
        streamflow, times, i_indices, j_indices = data_select.get_data_from_file(the_path)

        if not len(all_extremes):
            for i in range(streamflow.shape[1]):
                all_extremes.append([])


        for pos in range(streamflow.shape[1]):
            if high_flow:
                data1 = data_select.get_period_maxima(streamflow[:, pos], times,
                                start_date = start_date,
                                end_date = end_date,
                                start_month = start_month,
                                end_month = end_month,
                                event_duration = duration_days
                                )
            else:
                data1 = data_select.get_period_minima(streamflow[:, pos], times,
                                start_date = start_date,
                                end_date = end_date,
                                start_month = start_month,
                                end_month = end_month,
                                event_duration = duration_days
                                )
            all_extremes[pos].extend(list(data1.values()))

    #axes order: (time, position)
    all_extremes = np.array(all_extremes).transpose()

    if np.any(all_extremes is None):
        assert False, 'all_extremes = ' + str(all_extremes)

    #optimize
    print(all_extremes.shape)
    assert all_extremes.shape[1] == 547, 'all_extremes.shape[1] != 547'
    param_set = gevfit.optimize_stationary_for_period_and_all_cells_using_data(
                                        data = all_extremes,
                                        high_flow = high_flow
                                        )
    pickle.dump(param_set, open(param_file , 'wb'))
    return param_set
    pass

def fit_merged_for_current_and_future():
    """
    entry point to fit to all members at once
    """


    the_types = [True, False]

    high_start_month = 3
    high_end_month = 7
    high_duration = timedelta(days = 1)

    low_start_month = 1
    low_end_month = 5
    low_duration = timedelta(days = 15)

    type_to_period_start_month = {True : high_start_month, False : low_start_month}
    type_to_period_end_month = {True : high_end_month, False : low_end_month}

    type_to_event_duration = {True : high_duration, False : low_duration}


    data_folder = 'data/streamflows/hydrosheds_euler9'

    current_start_date = datetime(1970,1,1,0,0)
    current_end_date = datetime(1999,12, 31,0,0)

    future_start_date = datetime(2041,1,1,0,0)
    future_end_date = datetime(2070,12, 31,0,0)


    for hig_flow in the_types:
        pars_current = gev_fit_all_members(high_flow = hig_flow,
                        member_ids = members.current_ids,
                        data_folder = data_folder,
                        file_name_pattern = '{0}_discharge_1970_01_01_00_00.nc',
                        start_date = current_start_date,
                        end_date = current_end_date,
                        start_month = type_to_period_start_month[hig_flow] ,
                        end_month = type_to_period_end_month[hig_flow],
                        duration_days = type_to_event_duration[hig_flow])

        pars_future = gev_fit_all_members(high_flow = hig_flow,
                        member_ids = members.future_ids,
                        data_folder = data_folder,
                        file_name_pattern = '{0}_discharge_2041_01_01_00_00.nc',
                        start_date = future_start_date,
                        end_date = future_end_date,
                        start_month = type_to_period_start_month[hig_flow] ,
                        end_month = type_to_period_end_month[hig_flow],
                        duration_days = type_to_event_duration[hig_flow])



        assert len(pars_current) == 547
        assert len(pars_future) == 547

    pass




def main():
    fit_merged_for_current_and_future()



if __name__ == "__main__":
    main()
    print("Hello World")
