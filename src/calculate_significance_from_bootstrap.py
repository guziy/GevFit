__author__="huziy"
__date__ ="$3 fevr. 2011 09:35:44$"



import members
import gevfit

import pickle
import numpy as np
import data_select
import matplotlib.pyplot as plt
from shape.basin_boundaries import plot_basin_boundaries_from_shape
import matplotlib as mpl
from matplotlib.ticker import LinearLocator, MaxNLocator

import plot_utils

from math import sqrt
import pylab

import matplotlib_helpers.my_colormaps as mycolors



from map_parameters import polar_stereographic
basemap = polar_stereographic.basemap
xs = polar_stereographic.xs
ys = polar_stereographic.ys


inches_per_pt = 1.0 / 72.27               # Convert pt to inch
golden_mean = (sqrt(5) - 1.0) / 2.0       # Aesthetic ratio
fig_width = 600 * inches_per_pt          # width in inches
fig_height = fig_width * golden_mean      # height in inches
fig_size = [fig_width, 2.5 * fig_height]

font_size = 14

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


ret_level_getters = [
                     gevfit.get_high_ret_level_stationary,
                     gevfit.get_low_ret_level_stationary
                    ]

#level_type = 'high'/'low'
#returns list of lists of parameters for each point
#at each grid point: sigma, mu, ksi = pars
def get_pars_for_member_and_type(member_id, level_type = 'high'):
    file_name = 'gev_params_stationary_' + member_id + '_' + level_type
    return pickle.load(open(file_name))


##level_type = 'high'/'low'
def get_stdevs_for_member_and_type(member_id, level_type = 'high'):
    file_name = member_id + '_' + level_type + '_std_dev'
    return pickle.load(open(file_name))

def plot(data_1d, i_indices, j_indices, xs, ys,
         title = '', minmax = (None, None),
         color_map = mpl.cm.get_cmap('RdBu'),
         units = '',
         colorbar_orientation = 'vertical' , basemap = None,
         colorbar_tick_locator = LinearLocator(numticks = 6),
         colorbar_label_format = '%.1f', upper_limited = False, not_significant_mask = None):


    plt.title(title, {'fontsize': font_size})
    to_plot = np.ma.masked_all(xs.shape)
    sign_mask_2d = np.ma.masked_all(xs.shape)

    if not_significant_mask != None:
        the_zip = zip(i_indices, j_indices, data_1d, not_significant_mask)
        for i_index, j_index, the_data, significance in the_zip:
            to_plot[i_index, j_index] = the_data
            sign_mask_2d[i_index, j_index] = significance
        #plot not significant mask
        basemap.pcolormesh(xs, ys, sign_mask_2d.copy(), cmap = 'gray',
                     shading = 'flat', vmin = 0, vmax = 1 )
    else:
        the_zip = zip(i_indices, j_indices, data_1d)
        for i_index, j_index, the_data in the_zip:
            to_plot[i_index, j_index] = the_data





    basemap.pcolormesh(xs, ys, to_plot.copy(), cmap = color_map,
                    vmin = minmax[0], vmax = minmax[1], shading = 'flat', rasterized = False )


    plot_basin_boundaries_from_shape(basemap, plotter = plt, linewidth = 1, edge_color = 'k')
    basemap.drawcoastlines(linewidth = 0.5)
    
    plot_utils.draw_meridians_and_parallels(basemap, step_degrees = 30)



    ymin, ymax = plt.ylim()
    plt.ylim(ymin + 0.05 * (ymax - ymin) , ymax * 0.25)

    xmin, xmax = plt.xlim()
    plt.xlim(xmin + (xmax - xmin) * 0.55, 0.72*xmax)

    #plot colorbar
    cb = plt.colorbar(ticks = colorbar_tick_locator,
                      orientation = colorbar_orientation,
                      format = colorbar_label_format
                      )


    cb.ax.set_xlabel(units)
    if upper_limited:
        cl = cb.ax.get_yticklabels()
        labels = []
        for text in cl:
            labels.append(text.get_text())

        labels[-1] = '$\\geq$' + labels[-1]
        cb.ax.set_yticklabels(labels)






def calculate_and_plot(return_period = 10,
                       return_level_function = ret_level_getters[0]):

    if return_level_function == gevfit.get_high_ret_level_stationary:
        level_type = 'high'
    else:
        level_type = 'low'

    plt.clf()
    i_indices, j_indices = data_select.get_indices_from_file('data/streamflows/hydrosheds_euler9/aex_discharge_1970_01_01_00_00.nc')
    significance_counter = None
    plt.subplots_adjust(left = 0., hspace = 0.2, wspace = 0.2)

  

    for k, current_id in enumerate(members.current_ids):


        #get current return levels
        pars_list = get_pars_for_member_and_type(current_id, level_type)
        return_levels_current = np.zeros(len(pars_list))
        for pos, pars in enumerate(pars_list):
            return_levels_current[pos] = return_level_function(pars, return_period)
        stdevs_current = get_stdevs_for_member_and_type(current_id, level_type)[return_period]


        #get future return levels
        future_id = members.current2future[current_id]
        pars_list = get_pars_for_member_and_type(future_id, level_type)
        return_levels_future = np.zeros(len(pars_list))
        for pos, pars in enumerate(pars_list):
            return_levels_future[pos] = return_level_function(pars, return_period)
        stdevs_future = get_stdevs_for_member_and_type(future_id, level_type)[return_period]


        change = return_levels_future - return_levels_current
        if significance_counter == None:
            significance_counter = np.zeros( change.shape )


        print 'minmax(std_current)'
        print np.min(stdevs_current), np.max(stdevs_current)
        print 'minmax(std_future)'
        print np.min(stdevs_future), np.max(stdevs_future)

        print 'min max min abs(rl_current - rl_future)'
        the_delta = np.abs(return_levels_future - return_levels_current)
        print np.min(the_delta), np.max(the_delta), np.mean(the_delta)
        
        #stdev = -1 - stands for undefined value
        condition = np.logical_and(np.abs(change) > 1.96 * ( stdevs_current + stdevs_future ),
                                   (stdevs_current >= 0) & (stdevs_future >= 0)
                                   & (return_levels_current > 0)
                                )

        

       
        sign_index = np.where(condition)
        significance_counter[sign_index] += 1


        print len(sign_index[0])

        change[sign_index] /= return_levels_current[sign_index]
        change[sign_index] *= 100.0

        
        plt.subplot(3, 2, k + 1)
        if np.max(change[sign_index]) > 100:
            print np.max(change[sign_index])
            delta = 150
        else:
            delta = 50


        not_significant = np.ones(change.shape) * 0.75
        not_significant = np.ma.masked_where(condition, not_significant)


        if level_type == 'high':
            assert np.all(return_levels_current > 0)
            assert np.all(stdevs_current >= 0)
            assert np.all(stdevs_future >= 0)

        #temp change to condition
        change = np.ma.masked_where(np.logical_not(condition), change)
        print 'Plotting: current %s, future %s' % (current_id, future_id)
        plot(change , i_indices, j_indices, xs, ys,
                    title = '%s - %s' % (future_id, current_id),
                    color_map = mycolors.get_red_blue_colormap(ncolors = 16), units = '%',
                    basemap = basemap, minmax = (-delta, delta),
                    colorbar_label_format = '%d',
                    upper_limited = True, colorbar_tick_locator = LinearLocator(numticks = 9),
                    not_significant_mask = not_significant
                    )
      
        
    plt.subplot(3,2,6)

    significance_counter = np.ma.masked_where(significance_counter == 0, significance_counter)
    plot(significance_counter, i_indices, j_indices, xs, ys,
         title = 'Significance Count', minmax = (1,6),
         color_map = mycolors.get_sign_count_cmap(ncolors = 5), basemap = basemap,
         colorbar_tick_locator = MaxNLocator(nbins = 5),
         colorbar_label_format = '%d'
         )
    plt.savefig('%d_%s_change_rl.png' % (return_period, level_type), bbox_inches='tight')



def main():
    high_ret_periods = [10, 30, 50]
    low_ret_periods = [2, 5, 10]

    #temp for speed up
    high_ret_periods.pop()
    low_ret_periods.pop()

    for ret_period in high_ret_periods:
        calculate_and_plot(ret_period, gevfit.get_high_ret_level_stationary)
         

    for ret_period in low_ret_periods:
        calculate_and_plot(ret_period, gevfit.get_low_ret_level_stationary)
       

    pass


if __name__ == "__main__":
    main()
    print "Hello World"
