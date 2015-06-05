
__author__="huziy"
__date__ ="$31.12.2010 7:26:33$"

from multiprocessing import Pool
import os.path

import application_properties
import data_select
from datetime import datetime
from datetime import timedelta
import gevfit
import matplotlib as mpl
import members
import numpy as np
import numpy.random as random
import os
import pickle
application_properties.set_current_directory()


#np.seterr(all='raise', under='ignore')

def generate_indices(nvalues):
    return random.randint(0, nvalues, size = nvalues)

def generate_indices_restrict_data_to_member(nvalues, nvalues_per_member):
    result = []
    min_index = 0
    max_index = nvalues_per_member - 1
    n_members = nvalues / nvalues_per_member
    assert nvalues % nvalues_per_member == 0
    for i in range(n_members):
        result.extend(random.randint(min_index, high=max_index + 1, size=nvalues_per_member))
        min_index += nvalues_per_member
        max_index += nvalues_per_member
    return result

    pass

#
#return_periods - list of return periods
def function_for_each_process(args):
    sample_index, sampled_indices, extremes, \
    return_periods, high_flow, positions = args

    print('I work on sample %d' % sample_index)
    
    result = {}
    for the_period in return_periods:
        result[the_period] = []



    for pos in positions:
        pars = gevfit.optimize_stationary_for_period(extremes[sampled_indices, pos],
                                                    high_flow = high_flow)

        #calculate return levels for the current position and sample
        for the_period in return_periods:
            if high_flow:
                ret_level = gevfit.get_high_ret_level_stationary(pars, the_period)
            else:
                ret_level = gevfit.get_low_ret_level_stationary(pars, the_period)
            result[the_period].append(ret_level)
    return result



def apply_bootstrap_to_extremes(all_extremes, n_samples = 10,
                                out_file = "", process_pool = None,
                                return_periods = None, high_flow = True,
                                positions = None, restrict_indices_to_member = False,
                                n_values_per_member = -1
                                ):

    """
    calculate standard deviations using the bootstrap method
    and save to the binary file (dictionary, {ret_period => std dev array (1d)})
    all_extremes - values of extreme flow, axes: (time, position)

    """

    if os.path.isfile(out_file):
        print("{0} already exists, skipping ".format(out_file))
        return

    #prepare input data
    input = []
    for i in range(n_samples):
        if not restrict_indices_to_member:
            sampled_indices = generate_indices(all_extremes.shape[0])
        else:
            sampled_indices = generate_indices_restrict_data_to_member(all_extremes.shape[0], n_values_per_member)

        #input parameters
        input.append((i, sampled_indices, all_extremes,
                    return_periods, high_flow, positions))


    #returns list of maps which contain {return period : list of return levels
    #(the size of the list is npositions)}
    result = process_pool.map(function_for_each_process, input)

    assert len(result) == n_samples

    #standard deviations for all return periods
    all_std_devs = {}
    all_return_levels = {}

    for the_period in return_periods:
        all_return_levels[the_period] = []
        for the_dict in result:
            all_return_levels[the_period].append(the_dict[the_period])

    #calculate standard deviations
    for the_period in return_periods:
        sampled_ret_levels = np.array(all_return_levels[the_period])
        #all_std_devs[the_period] = - np.ones((len(positions),))


        #do not calculate the dispersion for the cases where some samples give
        #negative return levels (this is for the low flow return levels, high
        #flow return levels should always be positive)

        if not np.all(sampled_ret_levels >= 0):
            print("warning some resampled return levels were negative, assigning zeros")
            print("return period is %d " % the_period)
            sampled_ret_levels[sampled_ret_levels < 0] = 0.0

        all_std_devs[the_period] = np.std(sampled_ret_levels, axis = 0)

        print(np.array(all_return_levels[the_period]).shape)

    #sanity checks
    if high_flow:
        for the_period in return_periods:
            devs = all_std_devs[the_period]
            assert np.all(devs >= 0)
            assert np.all(np.array(all_return_levels) > 0)

    f = open(out_file, 'wb')
    pickle.dump(all_std_devs, f)
    f.close()



def apply_bootstrap_to_data(streamflow = None, times = None,
                    period_start_month = 1, period_end_month = 12,
                    start_date = datetime(1970,1,1,0,0),
                    end_date = datetime(1999,12,31,0,0),
                    event_duration_days = timedelta(days = 1),
                    n_samples = 1, high_flow = True,
                    return_periods = [], process_pool = None,
                    out_file = ''):


    """
    select extremes and perform bootstraping on streamflow
    the shape of the streamflow is (time, position)
    """
    
    all_extremes = []
    positions = [i for i in range(streamflow.shape[1])]

    for pos in positions:
        all_extremes.append([])

    for pos in positions:
        if high_flow:
            data = data_select.get_period_maxima(streamflow[:, pos], times,
                            start_date = start_date,
                            end_date = end_date,
                            start_month = period_start_month,
                            end_month = period_end_month,
                            event_duration = event_duration_days
                            )
        else:
            data = data_select.get_period_minima(streamflow[:, pos], times,
                            start_date = start_date,
                            end_date = end_date,
                            start_month = period_start_month,
                            end_month = period_end_month,
                            event_duration = event_duration_days
                            )
        all_extremes[pos].extend(list(data.values()))


    all_extremes = np.array(all_extremes).transpose()
    print('all_extremes.shape = ', all_extremes.shape)
    print('n_samples = ', n_samples)

    apply_bootstrap_to_extremes(all_extremes, n_samples = n_samples,
                                out_file = out_file,
                                return_periods = return_periods,
                                process_pool = process_pool,
                                positions = positions
                                )

    


#apply non-parametric bootstrap method in order to calculate
#uncertainties in return levels
#n_smaples - number of bootstrap samples
def apply_bootstrap(data_path = '',
                    member_name = 'aex',
                    period_start_month = 1, period_end_month = 12,
                    start_date = datetime(1970,1,1,0,0),
                    end_date = datetime(1999,12,31,0,0),
                    event_duration_days = timedelta(days = 1),
                    n_samples = 1, high_flow = True,
                    return_periods = [], process_pool = None
                    ):
    """
    Applying bootstrap to the given file at data_path
    """

    if high_flow:
        prefix = 'high'
    else:
        prefix = 'low'

    out_file = member_name + '_' + prefix + '_std_dev'
    if os.path.isfile(out_file):
        print('%s exists already ' % out_file)
        return


    #get streamflow data

    assert os.path.isfile(data_path)
    streamflow, times, xs, ys = data_select.get_data_from_file(data_path)
    
    apply_bootstrap_to_data(streamflow, times = times,
                        period_start_month = period_start_month,
                        period_end_month = period_end_month,
                        start_date = start_date,
                        end_date = end_date,
                        event_duration_days = event_duration_days,
                        n_samples = n_samples, high_flow = high_flow,
                        return_periods = return_periods,
                        process_pool = process_pool, out_file = out_file)





def plot(member = 'aex', low_return_periods = [], high_return_periods = []):

    for extreme in ['high', 'low']:
        data_path = member + '_' + extreme + '_std_dev'
        if not os.path.isfile(data_path):
            print('Warning: no data for (%s, %s) ' % (member, extreme))
            continue
        if not os.path.isfile(data_path):
            print('Warning can\'t find file %s for plotting' % data_path)
            continue

        if extreme == 'high':
            return_periods = high_return_periods
        else:
            return_periods = low_return_periods

        stds = pickle.load(open(data_path))
        for return_period in return_periods:
            imagefile = ( member + '_std_devs_%dyr_%s.png' ) % (return_period, extreme)
            gevfit.plot_data(data = stds[return_period],
                             imagefile = imagefile,
                             color_map = mpl.cm.get_cmap('RdBu',20))
    

import time
def main(n_samples = 20):
    low_return_periods = [2,5,10]

    high_return_periods = [10, 30, 50]


    #process pool
    n_processes = 12
    process_pool = Pool(processes = n_processes)


    #data_folder = 'data/streamflows/hydrosheds_euler9/'
    data_folder = "data/streamflows/narccap_ccsm-crcm"

    hi = "high"
    lo = "low"
    extreme_types = [hi, lo]

    high_start_month = 3
    high_end_month = 7
    high_duration = timedelta(days = 1)

    low_start_month = 1
    low_end_month = 5
    low_duration = timedelta(days = 15)

    #mappings from the type of extreme to the parameters
    type_to_period_start_month = {hi : high_start_month, lo : low_start_month}
    type_to_period_end_month = {hi : high_end_month, lo : low_end_month}
    type_to_event_duration = {hi : high_duration, lo : low_duration}
    type_to_return_periods = {hi: high_return_periods, lo : low_return_periods}


    current_start_date = datetime(1970,1,1,0,0)
    current_end_date = datetime(1999,12, 31,0,0)

    future_start_date = datetime(2041,1,1,0,0)
    future_end_date = datetime(2070,12, 31,0,0)


    #current_ids = members.current_ids
    #future_ids = members.future_ids
    current_ids = ["ccsm-crcm-current"]
    future_ids = ["ccsm-crcm-future"]


    plot_deviations = False
    for member in current_ids:
        data_path = '%s_discharge_1970_01_01_00_00.nc' % member
        data_path = os.path.join(data_folder, data_path)
        for extreme in extreme_types:
            apply_bootstrap(data_path = data_path,
                        member_name = member,
                        period_start_month = type_to_period_start_month[extreme],
                        period_end_month = type_to_period_end_month[extreme],
                        start_date = current_start_date,
                        end_date = current_end_date,
                        event_duration_days = type_to_event_duration[extreme],
                        n_samples = n_samples,
                        high_flow = (extreme == hi),
                        return_periods = type_to_return_periods[extreme],
                        process_pool = process_pool
                        )

        if plot_deviations:
            plot(member, low_return_periods = low_return_periods,
                         high_return_periods = high_return_periods)


    for member in future_ids:
        data_path = '%s_discharge_2041_01_01_00_00.nc' % member
        data_path = os.path.join(data_folder, data_path)

        for extreme in extreme_types:
            apply_bootstrap(data_path = data_path,
                        member_name = member,
                        period_start_month = type_to_period_start_month[extreme],
                        period_end_month = type_to_period_end_month[extreme],
                        start_date = future_start_date,
                        end_date = future_end_date,
                        event_duration_days = type_to_event_duration[extreme],
                        n_samples = n_samples,
                        high_flow = (extreme == hi),
                        return_periods = type_to_return_periods[extreme],
                        process_pool = process_pool
                        )

        if plot_deviations:
            plot(member, low_return_periods = low_return_periods,
                         high_return_periods = high_return_periods)
                         




if __name__ == "__main__":
    t0 = time.time()
    print(os.getcwd())
    main(n_samples = 1000)
    t1 = time.time()
    print('Execution time is %f seconds ' % (t1 - t0))
    print("Hello World")
