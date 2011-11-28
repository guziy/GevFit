__author__="huziy"
__date__ ="$Aug 29, 2011 4:25:10 PM$"

import members
from test_median import Test_median
import matplotlib.pyplot as plt
import application_properties
from datetime import datetime
from datetime import timedelta
import numpy as np

from scipy.stats import ttest_1samp
from mpl_toolkits.basemap import Basemap

import pickle
import os
from map_parameters import polar_stereographic



def get_test_median_objects():

    test_median_data_path = 'test_median.bin'
    if os.path.isfile(test_median_data_path):
        test_median_current, test_median_future = pickle.load(open(test_median_data_path))
    else:
        paths = [
            'data/streamflows/hydrosheds_euler9/aet_discharge_1970_01_01_00_00.nc',
            'data/streamflows/hydrosheds_euler9/aeu_discharge_2041_01_01_00_00.nc'
        ]



        path_to_start_date = {paths[0]:datetime(1970,1,1), paths[1]:datetime(2041,1,1)}

        start_to_end = {datetime(1970,1,1): datetime(1999,12, 31), datetime(2041,1,1) : datetime(2070,12, 31)}


        type_to_startmonth = {'low':1,'high':3}
        type_to_end_month = {'low':5,'high':7}
        type_to_duration = {'low':timedelta(days = 15),'high':timedelta(days = 1)}

        #low_return_periods = [2, 5]

        the_type = 'high'
        the_period = 2 #years

        #current data
        test_median_current = Test_median(data_path = paths[0])
        start_date = path_to_start_date[paths[0]]
        test_median_current.start_date = start_date
        test_median_current.end_date = start_to_end[start_date]
        test_median_current.start_month = type_to_startmonth[the_type]
        test_median_current.end_month = type_to_end_month[the_type]
        test_median_current.return_period_years = the_period
        test_median_current.high_flow = True
        test_median_current.event_duration = type_to_duration[the_type]
        test_median_current.select_data_and_calculate()

        #future data
        test_median_future = Test_median(data_path = paths[1])
        start_date = path_to_start_date[paths[1]]
        test_median_future.start_date = start_date
        test_median_future.end_date = start_to_end[start_date]
        test_median_future.start_month = type_to_startmonth[the_type]
        test_median_future.end_month = type_to_end_month[the_type]
        test_median_future.return_period_years = the_period
        test_median_future.high_flow = True
        test_median_future.event_duration = type_to_duration[the_type]
        test_median_future.select_data_and_calculate()

        pickle.dump( [test_median_current, test_median_future], open(test_median_data_path , 'w'))

    return test_median_current, test_median_future




def main():

    test_median_current, test_median_future = get_test_median_objects()
    c_median = test_median_current.median_field
    f_median = test_median_future.median_field

    c_rl = test_median_current.ret_level_2yr
    f_rl = test_median_future.ret_level_2yr

    median_changes = (f_median - c_median) / c_median * 100.0
    rl_changes = (f_rl - c_rl) / c_rl * 100.0

    #compare median and 2 year return level
    plt.subplot(1,2,1)
    plt.title('values in ${\\rm m^3/s}$,\n current climate')
    plt.scatter(c_median, c_rl)
    x0 = min( np.min(c_median), np.min(c_rl) )
    x1 = max( np.max(c_median), np.max(c_rl) )
    plt.plot([x0, x1], [x0, x1], color = 'k')
    plt.xlabel('median')
    plt.ylabel('2 year return level')

    plt.subplot(1,2,2)
    plt.title('values in ${\\rm m^3/s}$, \n future climate')
    plt.scatter(f_median, f_rl)
    x0 = min( np.min(f_median), np.min(f_rl) )
    x1 = max( np.max(f_median), np.max(f_rl) )
    plt.plot([x0, x1], [x0, x1], color = 'k')
    plt.xlabel('median')
    plt.ylabel('2 year return level')



    #compare the changes in median and 2 year return level
    plt.figure()
    plt.subplots_adjust(hspace = 0.3)


    plt.subplot(1,2,1)
    plt.scatter(median_changes, rl_changes)
    x0 = min( np.min(median_changes), np.min(rl_changes) )
    x1 = max( np.max(median_changes), np.max(rl_changes) )
    plt.plot([x0, x1], [x0, x1], color = 'k')
    plt.xlabel('median')
    plt.ylabel('2 year return level')
    plt.title('changes')
    
    plt.subplot(1,2,2)
    pvalues = np.zeros((median_changes.shape[0],))

    for i in xrange(len(pvalues)):
        t, pvalues[i] = ttest_1samp( np.array(test_median_current.data_extremes[i]) -
                                     np.array(test_median_future.data_extremes[i]), 0)


    # @type test_median_current Test_median
    print len(pvalues), test_median_current.longitudes.shape

    b = Basemap(resolution = 'i')
    i_indices, j_indices = test_median_current.get_indices_in_2d_grid()
    to_plot = np.ma.masked_all(polar_stereographic.lons.shape)
    for i, j, pvalue in zip(i_indices, j_indices, pvalues):
        to_plot[i,j] = pvalue

    b.pcolormesh(polar_stereographic.lons, polar_stereographic.lats, to_plot)
    b.drawcoastlines()
    plt.colorbar(shrink = 0.5)


    plt.title('p-value of student \n ttest applied to the extree data')
    plt.show()



if __name__ == "__main__":
    application_properties.set_current_directory()
    main()
    print "Hello World"
