from datetime import timedelta, datetime
import itertools
import pylab
import plot_utils

__author__ = 'huziy'

import data_select
import os
import numpy as np
from scipy.stats.mstats import kruskalwallis
import calculate_significance_from_bootstrap as csfb
import members
import matplotlib.pyplot as plt
import matplotlib as mpl

from map_parameters import polar_stereographic

import application_properties
import pickle




def get_extremes_list(data_path = "", member_ids = None, high_flow = True,
                        start_date = None, end_date = None,
                        event_duration = timedelta(days = 1),
                        period_start_month = 1, period_end_month = 12
                        ):
    """
    returns list of 2d arrays of extremes, the 2d arrays have the  shape = (time, cell_index)
    """
    file_paths = []
    for the_name in os.listdir(data_path):
        prefix = the_name.split('_')[0]
        if prefix in member_ids:
            file_paths += [os.path.join(data_path, the_name)]


    #merge extreme data
    all_extremes = []
    i_indices = None
    j_indices = None
    for the_path in file_paths:
        streamflow, times, i_indices, j_indices = data_select.get_data_from_file(the_path)


        domain_extremes = [[] for pos in range(len(i_indices))]

        for pos, point_extrems in enumerate(domain_extremes):
            if high_flow:
                extremes = data_select.get_period_maxima(streamflows=streamflow[:, pos], times = times,
                                                               start_date = start_date, end_date = end_date,
                                                               event_duration = event_duration,
                                                               start_month = period_start_month,
                                                               end_month = period_end_month
                                                               )
            else:
                extremes = data_select.get_period_minima(streamflows=streamflow[:, pos], times = times,
                                                           start_date = start_date, end_date = end_date,
                                                           event_duration = event_duration,
                                                           start_month = period_start_month,
                                                           end_month = period_end_month
                                                           )
            point_extrems.extend(list(extremes.values()))

        all_extremes.append(np.transpose( np.array(domain_extremes) ))


    return all_extremes, i_indices, j_indices


def kw_test(data):
    """
    data = [data2d_1,..,data2d_n]
    """
    n_pos = data[0].shape[1]
    p_values = np.zeros((n_pos,))

    for pos in range(n_pos):
        samples = [
            data2d[:, pos] for data2d in data
        ]
        h, p_values[pos] = kruskalwallis(*samples)

    return p_values



class KeyDataObject:
    def __init__(self, the_type = "high", time_window = "current", p_values = None):
        self.type = the_type
        self.time_window = time_window
        self.p_values = p_values


def get_key_data_list():

    dump_file = "kw_test_data.bin"
    if os.path.isfile(dump_file):
        return pickle.load(open(dump_file))


    hi = "high"
    lo = "low"

    the_types = [hi, lo]

    high_start_month = 3
    high_end_month = 7
    high_duration = timedelta(days = 1)

    low_start_month = 1
    low_end_month = 5
    low_duration = timedelta(days = 15)

    type_to_period_start_month = {hi : high_start_month, lo : low_start_month}
    type_to_period_end_month = {hi : high_end_month, lo : low_end_month}

    type_to_event_duration = {hi : high_duration, lo : low_duration}


    data_folder = 'data/streamflows/hydrosheds_euler9'

    current_start_date = datetime(1970,1,1,0,0)
    current_end_date = datetime(1999,12, 31,0,0)

    future_start_date = datetime(2041,1,1,0,0)
    future_end_date = datetime(2070,12, 31,0,0)

    fc = "future"
    cc = "current"
    time_windows = [cc, fc]
    tw_to_ids_list = {
        cc: members.current_ids,
        fc: members.future_ids
    }
    key_data_list = []
    for tw in time_windows:
        mlist = tw_to_ids_list[tw]

        start_date = current_start_date if tw  == cc else future_start_date
        end_date = current_end_date if tw  == cc else future_end_date

        for the_type in the_types:
            start_month = type_to_period_start_month[the_type]
            end_month = type_to_period_end_month[the_type]
            the_duration = type_to_event_duration[the_type]
            extremes, i_list, j_list = get_extremes_list(data_path=data_folder,
                                         member_ids = mlist, high_flow = (the_type == hi),
                                         event_duration = the_duration,
                                         start_date=start_date, end_date=end_date,
                                         period_start_month=start_month,
                                         period_end_month=end_month
            )

            p_values = kw_test(extremes)
            kv = KeyDataObject(the_type=the_type, time_window=tw, p_values=p_values)
            key_data_list.append(kv)
    pickle.dump(key_data_list, open(dump_file, "w"))
    return key_data_list


def kw_test_for_means(current_climate = True, data_folder = 'data/streamflows/hydrosheds_euler9', months = list(range(1,13))):
    """
    returns p-values resulting from kruskal - wallis test on annual means
    """

    the_ids = members.all_current if current_climate else members.all_future

    file_paths = []
    for the_file in os.listdir(data_folder):
        if the_file.split("_")[0] in the_ids:
            file_paths.append(os.path.join(data_folder, the_file))

    real_means = []
    for the_path in file_paths:
        streamflow, times, i_indices, j_indices = data_select.get_data_from_file(the_path)

        #for each year and for each gridcell get mean value for the period
        means_dict = data_select.get_means_over_months_for_each_year(times, streamflow, months = months)

        means_sorted_in_time = [x[1] for x in sorted(list(means_dict.items()), key=lambda x: x[0])]
        data_matrix = np.array(means_sorted_in_time)
        real_means.append(data_matrix) #save modelled means
        #print "data_matrix.shape = ", data_matrix.shape

    n_positions = real_means[0].shape[1]
    p_values = np.zeros((n_positions,))
    for pos in range(n_positions):
        samples = [
            data2d[:, pos] for data2d in real_means
        ]

        #x = list(samples)
        #print len(x), x[0].shape


        h, p_values[pos] = kruskalwallis(*samples)
    return p_values

    pass



def main():
    plot_utils.apply_plot_params(width_pt=None, font_size=9, aspect_ratio=2.5)
    gs = mpl.gridspec.GridSpec(3,2)
    key_data_list = get_key_data_list()
    i_list, j_list = data_select.get_indices_from_file()
    subplot_count = 0

    for the_type in ["high", "low"]:
        for time_window in ["current", "future"]:
            selected_data = None

            for data in key_data_list:
                if data.time_window == time_window and data.type == the_type:
                    selected_data = data
                    break

            row = subplot_count // 2
            col = subplot_count % 2
            plt.subplot(gs[row, col])
            csfb.plot(selected_data.p_values, i_list, j_list,
                      polar_stereographic.xs, polar_stereographic.ys,
                      units = "", basemap = polar_stereographic.basemap,
                      minmax = (0, 0.25), title = "", # "{0} climate, {1} flow".format(selected_data.time_window, selected_data.type),
                      colorbar_label_format="%.2f", color_map = mpl.cm.get_cmap("jet", 5), upper_limited=True
            )
            subplot_count += 1

    #TODO:add 2 subplots for mean values
    pc = kw_test_for_means()
    plt.subplot(gs[2,0])
    csfb.plot(pc, i_list, j_list,
              polar_stereographic.xs, polar_stereographic.ys,
              units = "", basemap = polar_stereographic.basemap,
              minmax = (0, 0.25), #title = "{0} climate, {1} flow".format("current", "mean"),
              colorbar_label_format="%.2f", color_map = mpl.cm.get_cmap("jet", 5), upper_limited=True
    )

    pf = kw_test_for_means(current_climate=False)
    plt.subplot(gs[2, 1])
    csfb.plot(pf, i_list, j_list,
              polar_stereographic.xs, polar_stereographic.ys,
              units = "", basemap = polar_stereographic.basemap,
              minmax = (0, 0.25), #title = "{0} climate, {1} flow".format("future", "mean"),
              colorbar_label_format="%.2f", color_map = mpl.cm.get_cmap("jet", 5), upper_limited=True
    )


    plt.tight_layout()
    plt.savefig("p_values_kruskalwallis.png")

def test():
    application_properties.set_current_directory()
    main()


if __name__ == "__main__":
    test()