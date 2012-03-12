from matplotlib import gridspec
import plot_utils

__author__="huziy"
__date__ ="$Sep 17, 2011 6:44:52 PM$"


from datetime import datetime
from datetime import timedelta
import data_select
import bootstrap
import application_properties

import members

from multiprocessing import Pool
import os

import pickle
import numpy as np
import gevfit

from map_parameters import polar_stereographic

import calculate_significance_from_bootstrap as csfb
import matplotlib_helpers.my_colormaps as mycolors
from matplotlib.ticker import LinearLocator

import matplotlib.pyplot as plt

def apply_bootstrap_to_all_members_merged(file_paths = None,
                                high_flow = True,
                                n_samples = 10, out_file = '',
                                process_pool = None,
                                start_date = None,
                                end_date = None,
                                start_month = None,
                                end_month = None,
                                duration_days = None,
                                return_periods = None
                                ):
    """
    duration_days - timedelta object
    """

    if os.path.isfile(out_file):
        print "{0} already exists, skipping ".format(out_file)
        return



    #select data
    all_extremes = []
    streamflow = None
    for the_path in file_paths:
        print the_path
        streamflow, times, i_indices, j_indices = data_select.get_data_from_file(the_path)

        if not len(all_extremes):
            all_extremes = [[] for i in xrange(streamflow.shape[1])]

        for pos in xrange(streamflow.shape[1]):
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
            all_extremes[pos].extend(data1.values())

    #axes order: (time, position)
    all_extremes = np.array(all_extremes).transpose()
    bootstrap.apply_bootstrap_to_extremes(all_extremes,
                                        n_samples = n_samples,
                                        out_file = out_file,
                                        process_pool = process_pool,
                                        return_periods = return_periods,
                                        positions = xrange(streamflow.shape[1]),
                                        high_flow = high_flow,
                                        restrict_indices_to_member=True,
                                        n_values_per_member= all_extremes.shape[0] / len(file_paths)
                                        )
    print "n_indices per member = ", all_extremes.shape[0] / len(file_paths)
    pass



def get_high_return_periods():
    return [10, 30]

def get_low_return_periods():
    return [2, 5]

def main( n_samples = 10 ):
    data_folder = 'data/streamflows/hydrosheds_euler9'

    hi = "high"
    lo = "low"
    extreme_types = [ hi, lo]

    high_start_month = 3
    high_end_month = 7
    high_duration = timedelta(days = 1)
    high_return_periods = get_high_return_periods()

    low_start_month = 1
    low_end_month = 5
    low_duration = timedelta(days = 15)
    low_return_periods = get_low_return_periods()

    type_to_period_start_month = {hi : high_start_month, lo : low_start_month}
    type_to_period_end_month = {hi : high_end_month, lo : low_end_month}
    type_to_event_duration = {hi : high_duration, lo : low_duration}
    type_to_return_periods = {hi : high_return_periods, lo : low_return_periods}


    current_start_date = datetime(1970,1,1,0,0)
    current_end_date = datetime(1999,12, 31,0,0)

    future_start_date = datetime(2041,1,1,0,0)
    future_end_date = datetime(2070,12, 31,0,0)

    #get paths of the data files to be merged for bootstrap
    current_paths = []
    future_paths = []
    for fName in os.listdir(data_folder):
        the_id = fName.split('_')[0]

        if the_id in members.current_ids:
            current_paths.append(os.path.join(data_folder, fName))
            continue

        if the_id in members.future_ids:
            future_paths.append(os.path.join(data_folder, fName))



    #process pool
    #n_processes = int(n_samples ** 0.6) if n_samples > 100 else n_samples
    process_pool = Pool(processes = 20)

    for extreme in extreme_types:
        start_month = type_to_period_start_month[extreme]
        end_month = type_to_period_end_month[extreme]
        event_duration = type_to_event_duration[extreme]
        return_periods = type_to_return_periods[extreme]

        #current
        apply_bootstrap_to_all_members_merged(
                    file_paths = current_paths,
                    high_flow = (extreme == hi),
                    n_samples = n_samples,
                    out_file = extreme + "_std_dev_current",
                    process_pool = process_pool,
                    start_date = current_start_date,
                    end_date = current_end_date,
                    start_month = start_month,
                    end_month = end_month,
                    duration_days = event_duration,
                    return_periods = return_periods
                    )
        #future
        apply_bootstrap_to_all_members_merged(
                    file_paths = future_paths,
                    high_flow = (extreme == hi),
                    n_samples = n_samples,
                    out_file = extreme + "_std_dev_future",
                    process_pool = process_pool,
                    start_date = future_start_date,
                    end_date = future_end_date,
                    start_month = start_month,
                    end_month = end_month,
                    duration_days = event_duration,
                    return_periods = return_periods
                    )
    pass


def plot_results():

    cc = "current"
    fc = "future"
    climate_types = [cc, fc]


    folder_path = 'data/streamflows/hydrosheds_euler9/'
    file_path = os.path.join(folder_path, "aex_discharge_1970_01_01_00_00.nc")
    i_indices, j_indices = data_select.get_indices_from_file(file_path)
    xs = polar_stereographic.xs
    ys = polar_stereographic.ys
    basemap = polar_stereographic.basemap
    
    hi = "high"
    lo = "low"
    extreme_types = [hi, lo]

    base_std_name = "{0}_std_dev_{1}"
    base_par_name = {
        cc : "_".join(["{0}"] + members.current_ids),
        fc : "_".join(["{0}"] + members.future_ids)
    }

    extreme_to_return_periods = {
        hi : get_high_return_periods(),
        lo : get_low_return_periods()
    }


    sig_coefs = [1.96
    #    , 1.645
    ]
    sig_levels = ["95 %"
    #    , "90 %"
    ]
    gs = gridspec.GridSpec(3,2)

    for extreme in extreme_types:
        pars_path_current = base_par_name[cc].format(extreme)
        std_path_current = base_std_name.format(extreme, cc)

        pars_path_future = base_par_name[fc].format(extreme)
        std_path_future = base_std_name.format(extreme, fc)

        pars_current = pickle.load(open(pars_path_current))
        stds_current = pickle.load(open(std_path_current))

        pars_future = pickle.load(open(pars_path_future))
        stds_future = pickle.load(open(std_path_future))

        return_periods = extreme_to_return_periods[extreme]
        #plot_utils.apply_plot_params(font_size=15, width_pt=900, aspect_ratio=2.5)
        #plt.figure()


        delta = 50 if extreme == hi else 100
        for row, ret_period in enumerate( return_periods ):
            #calculate changes in return levels
            func = lambda x: gevfit.get_high_ret_level_stationary(x, ret_period)
            rl_c = map(func, pars_current)
            rl_f = map(func, pars_future)
            
            rl_c = np.array(rl_c)
            rl_f = np.array(rl_f)

            std_c = stds_current[ret_period]
            std_f = stds_future[ret_period]

            in_odf = (std_c > 0) & (std_f > 0) & (rl_c > 0) & (rl_f >= 0)
            change = np.ma.masked_all(rl_c.shape)
            change[in_odf] = (rl_f[in_odf] - rl_c[in_odf]) / rl_c[in_odf] * 100.0


            min_change = np.min((rl_f - rl_c) / rl_c * 100.0)
            if min_change >= 0:
               low_limit = 0
            elif min_change > -10:
               low_limit = -10
            else:
                low_limit = np.floor(min_change / 10.0) * 10

            print "min change = {0}, low limit = {1}".format(min_change, low_limit)


            for sig_coef, sig_name in zip(sig_coefs, sig_levels):
                significance = np.ma.masked_all(rl_c.shape)
                sig_cond = (sig_coef * (std_c + std_f) < np.abs(rl_f - rl_c)) & in_odf
                significance[~sig_cond] = 0#fill with gray non signifcant areas
                change1 = np.ma.masked_where(~sig_cond, change)
                if extreme == hi:
                    plt.subplot(gs[row, 0])
                else:
                    plt.subplot(gs[row, 1])

                csfb.plot(change1 , i_indices, j_indices, xs, ys,
                        title = 'T = {0}-year'.format( ret_period ),
                        color_map = mycolors.get_red_blue_colormap(ncolors = 10),
                        units = '%',
                        basemap = basemap, minmax = (-delta, delta),
                        colorbar_label_format = '%d',
                        upper_limited = True,
                        colorbar_tick_locator = LinearLocator(numticks = 11),
                        not_significant_mask = significance,
                        show_colorbar = True, impose_lower_limit=low_limit
                        )
                #subplot_count += 1
    plt.tight_layout()
    plt.savefig("rl_of_merged_change.png")

    pass

if __name__ == "__main__":
    application_properties.set_current_directory()
    plot_utils.apply_plot_params(width_pt=None, font_size=9)
    calculate = False
    if calculate:
        main(n_samples = 1000)
    else:
        plot_results()
    print "Hello World"
