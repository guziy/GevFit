import os
import data_select
import gevfit
from datetime import datetime, timedelta

import numpy as np

from map_parameters import polar_stereographic
import matplotlib.pyplot as plt
from shape.basin_boundaries import plot_basin_boundaries_from_shape

import application_properties
application_properties.set_current_directory()

class Test_median():
    def __init__(self, data_path = ''):
        data, times, i_indices, j_indices = data_select.get_data_from_file(data_path)
        self._id, rest = os.path.basename(data_path).split('_', 1)
        self._data = data
        self._times = times
        self._i_indices = i_indices
        self._j_indices = j_indices

        self.data_extremes = None
        self.return_period_years = 2
        self.high_flow = True #low_flows are calculated if False        

        self.start_date = None
        self.end_date = None
        self.start_month = 1
        self.end_month = 12

        self.event_duration = timedelta(days = 1)

        self.median_field = None
        self.ret_level_2yr = None



    def select_data_and_calculate(self):
        if self.high_flow:
            self.data_extremes = data_select.get_list_of_annual_maximums_for_domain(self._data, self._times,
                                        start_date = self.start_date, end_date = self.end_date,
                                        start_month = self.start_month, end_month = self.end_month,
                                        event_duration = self.event_duration)
        else:
            self.data_extremes = data_select.get_list_of_annual_minimums_for_domain(self._data, self._times,
                                        start_date = self.start_date, end_date = self.end_date,
                                        start_month = self.start_month, end_month = self.end_month,
                                        event_duration = self.event_duration)


        the_type = 'high' if self.high_flow else 'low'
        save_extremes_to_txt_file('{0}_{1}_values.txt'.format(self._id, the_type), self.data_extremes, self._i_indices, self._j_indices)
        
        self._calculate_median_field()
        print 'calculated median field'
        self._calculate_return_level_field()
        pass

    def _calculate_median_field(self):
        assert self.start_date != None and self.end_date != None, 'start_date and end_date fields should be set.'
        self.median_field = []
        #cycle through the points of the domain
        for pos in range(len(self._i_indices)):
            values = self.data_extremes[pos]
            if np.median(values) > 10000:
                print values
            self.median_field.append(np.median(values))
        pass

    def _calculate_return_level_field(self):
        if self.high_flow:
            field = gevfit.get_high_levels_for_id(self._id, return_period = self.return_period_years)
        else:
            field = gevfit.get_low_levels_for_id(self._id, return_period = self.return_period_years)



        the_type = 'high' if self.high_flow else 'low'

        save_ret_levels_to_txt('{0}_{1}yr_{2}_ret_level.txt'.format(self._id, self.return_period_years, the_type),
                                    field, self._i_indices, self._j_indices)

        save_pars_to_txt_file('{0}_{1}_params.txt'.format(self._id, the_type),
                            gevfit.get_gevd_params_for_id_and_type(self._id, self.high_flow) ,self._i_indices, self._j_indices)



        self.ret_level_2yr = []
        for k in range(len(self._i_indices)):
            self.ret_level_2yr.append(field[k])


    def plot(self):
        basemap = polar_stereographic.basemap

        plt.subplot(3,1,1)
        gevfit.plot_data(self.ret_level_2yr, imagefile = None,
                            units = 'm**3/s', minmax = (0, None),
                            i_list = self._i_indices, j_list = self._j_indices)
        plt.title('2 year return levels')

        plt.subplot(3,1,2)
        gevfit.plot_data(self.median_field, imagefile = None,
                            units = 'm**3/s', minmax = (0, None),
                            i_list = self._i_indices, j_list = self._j_indices)
        plt.title('median')
        
        
        plt.subplot(3,1,3)
        gevfit.plot_data(np.array(self.median_field) - np.array(self.ret_level_2yr), imagefile = None,
                            units = 'm**3/s', minmax = (None, None),
                            i_list = self._i_indices, j_list = self._j_indices)

        plt.title('difference: median - 2 year return level')
        plot_basin_boundaries_from_shape(basemap, plotter = plt, linewidth = 0.1)
        basemap.drawcoastlines()

        plt.savefig('median_test.png')
        pass


def save_pars_to_txt_file(filename, pars_dict, i_indices, j_indices):
    f = open(filename, 'w')
    f.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format('i', 'j', 'sigma', 'mu', 'ksi'))
    for pos in range(len(i_indices)):
        f.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(i_indices[pos], j_indices[pos],
                            pars_dict[pos][0], pars_dict[pos][1], pars_dict[pos][2]))
    f.close()
    pass

def save_extremes_to_txt_file(filename, data, i_indices, j_indices):
    f = open(filename, 'w')
    for i, j, the_value_list in zip(i_indices, j_indices, data):
        f.write('i = {0},\t j = {1}\n'.format(i, j))
        for value in the_value_list:
            f.write('{0}\n'.format(value))
    f.close()
    pass

def save_ret_levels_to_txt(filename, data, i_indices, j_indices):
    f = open(filename, 'w')
    for i, j, the_value in zip(i_indices, j_indices, data):
        f.write('i = {0},\t j = {1}, \t {2} \n'.format(i, j, the_value))
    f.close()



def main():


    paths = [
        'data/streamflows/hydrosheds_euler3/aet_discharge_1970_01_01_00_00.nc',
        'data/streamflows/hydrosheds_euler3/aeu_discharge_2041_01_01_00_00.nc'
    ]


    path_to_start_date = {paths[0]:datetime(1970,1,1), paths[1]:datetime(2041,1,1)}

    start_to_end = {datetime(1970,1,1): datetime(1999,12, 31), datetime(2041,1,1) : datetime(2070,12, 31)}

    the_types = ['low', 'high']

    type_to_startmonth = {'low':3,'high':4}
    type_to_end_month = {'low':4,'high':6}
    type_to_duration = {'low':timedelta(days = 15),'high':timedelta(days = 1)}

    low_return_periods = [2, 5]
    high_return_periods = [30, 50, 100]

    the_type = 'high'
    for path in paths:
        for the_period in high_return_periods:
            test_median = Test_median(data_path = path)
            start_date = path_to_start_date[path]
            test_median.start_date = start_date
            test_median.end_date = start_to_end[start_date]
            test_median.start_month = type_to_startmonth[the_type]
            test_median.end_month = type_to_end_month[the_type]
            test_median.return_period_years = the_period
            test_median.high_flow = True
            test_median.event_duration = type_to_duration[the_type]
            print 'selecting data and calculating'

            print 'path={0}, the_type = {1}, the_period = {2}'.format(path, the_type, the_period)
            test_median.select_data_and_calculate()

    the_type = 'low'
    for path in paths:
        for the_period in low_return_periods:
            test_median = Test_median(data_path = path)
            start_date = path_to_start_date[path]
            test_median.start_date = start_date
            test_median.end_date = start_to_end[start_date]
            test_median.start_month = type_to_startmonth[the_type]
            test_median.end_month = type_to_end_month[the_type]
            test_median.return_period_years = the_period
            test_median.high_flow = False
            test_median.event_duration = type_to_duration[the_type]
            print 'selecting data and calculating'

            print 'path={0}, the_type = {1}, the_period = {2}'.format(path, the_type, the_period)
            test_median.select_data_and_calculate()





#    data_path = 'data/streamflows/hydrosheds_euler2_1/aet_discharge_1970_01_01_00_00.nc'
#    test_median = Test_median(data_path = data_path)
#    test_median.start_date = datetime(1970,1,1)
#    test_median.end_date = datetime(1999,12, 31)
#    test_median.start_month = 4
#    test_median.end_month = 6
#    test_median.return_period_years = 2
#    test_median.high_flow = True
#    test_median.select_data_and_calculate()
#    test_median.plot()

#    data_path = 'data/streamflows/hydrosheds_euler2/aeu_discharge_2041_01_01_00_00.nc'
#    test_median = Test_median(data_path = data_path)
#    test_median.start_date = datetime(2041,1,1)
#    test_median.end_date = datetime(2070,12, 31)
#    test_median.start_month = 4
#    test_median.end_month = 6
#    test_median.return_period_years = 30
#    test_median.select_data_and_calculate()
#    test_median.plot()
#



if __name__ == '__main__':
    main()
    pass

