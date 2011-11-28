from datetime import timedelta, datetime
import pylab

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


def apply_plot_params(font_size = 16):
    """

    """
    inches_per_pt = 1.0 / 72.27               # Convert pt to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0       # Aesthetic ratio
    fig_width = 900 * inches_per_pt          # width in inches
    fig_height = fig_width * golden_mean      # height in inches
    fig_size = [fig_width, 1.5 * fig_height]


    params = {
            'axes.labelsize': font_size,
            'font.size':font_size,
            'text.fontsize': font_size,
            'legend.fontsize': font_size,
            'xtick.labelsize': font_size,
            'ytick.labelsize': font_size,
            'figure.figsize': fig_size
            }

    pylab.rcParams.update(params)



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


        domain_extremes = [[] for pos in xrange(len(i_indices))]

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
            point_extrems.extend(extremes.values())

        all_extremes.append(np.transpose( np.array(domain_extremes) ))


    return all_extremes, i_indices, j_indices


def kw_test(data):
    """
    data = [data2d_1,..,data2d_n]
    """
    n_pos = data[0].shape[1]
    p_values = np.zeros((n_pos,))

    for pos in xrange(n_pos):
        samples = (
            data2d[:, pos] for data2d in data
        )
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



def main():
    apply_plot_params()
    key_data_list = get_key_data_list()
    i_list, j_list = data_select.get_indices_from_file()
    subplot_count = 1
    for data in key_data_list:
        plt.subplot(2,2, subplot_count)
        csfb.plot(data.p_values, i_list, j_list,
                  polar_stereographic.xs, polar_stereographic.ys,
                  units = "", basemap = polar_stereographic.basemap,
                  minmax = (0, 0.2), title = "{0} climate, {1} flow".format(data.time_window, data.type),
                  colorbar_label_format="%.2f", color_map = mpl.cm.get_cmap("jet", 5), upper_limited=True
        )
        subplot_count += 1

    plt.savefig("p_values_kruskalwallis.pdf", bbox_inches = "tight")

def test():
    application_properties.set_current_directory()
    main()


if __name__ == "__main__":
    test()