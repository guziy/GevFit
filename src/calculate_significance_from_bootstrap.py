__author__="huziy"
__date__ ="$3 fevr. 2011 09:35:44$"



from data.query_object import QueryObject
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

from mpl_toolkits.basemap import Basemap


from datetime import datetime
from datetime import timedelta

from map_parameters import polar_stereographic
basemap = polar_stereographic.basemap
xs = polar_stereographic.xs
ys = polar_stereographic.ys



import save_to_file_rls_and_sign as txt_saver

from matplotlib.colors import ListedColormap


inches_per_pt = 1.0 / 72.27               # Convert pt to inch
golden_mean = (sqrt(5) - 1.0) / 2.0       # Aesthetic ratio
fig_width = 900 * inches_per_pt          # width in inches
fig_height = fig_width * golden_mean      # height in inches
fig_size = [fig_width, 2.5 * fig_height]

font_size = 13

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
         title = '', label = "", minmax = (None, None),
         color_map = mpl.cm.get_cmap('RdBu'),
         units = '',
         colorbar_orientation = 'vertical' , basemap = None,
         colorbar_tick_locator = LinearLocator(numticks = 6),
         colorbar_label_format = '%.1f', upper_limited = False,
         not_significant_mask = None, show_colorbar = True):


    plt.title(title, {'fontsize': font_size})
    to_plot = np.ma.masked_all(xs.shape)
    sign_mask_2d = np.ma.masked_all(xs.shape)

    if not_significant_mask is not None:
        the_zip = zip(i_indices, j_indices, data_1d, not_significant_mask)
        for i_index, j_index, the_data, significance in the_zip:
            to_plot[i_index, j_index] = the_data
            sign_mask_2d[i_index, j_index] = significance
        #plot not significant mask
        basemap.pcolormesh(xs, ys, sign_mask_2d.copy(), cmap = mpl.cm.get_cmap('gist_gray', 3),
                     shading = 'flat', vmin = -1, vmax = 1, zorder = 2 )
    else:
        the_zip = zip(i_indices, j_indices, data_1d)
        for i_index, j_index, the_data in the_zip:
            to_plot[i_index, j_index] = the_data




    plot_axes = plt.gca()
    image = basemap.pcolormesh(xs, ys, to_plot.copy(), cmap = color_map,
                    vmin = minmax[0], vmax = minmax[1], shading = 'flat', rasterized = False )


    plot_basin_boundaries_from_shape(basemap, plotter = plt, linewidth = 1, edge_color = 'k')
    basemap.drawcoastlines(linewidth = 0.5)
    
    plot_utils.draw_meridians_and_parallels(basemap, step_degrees = 30)



    x_min, x_max, y_min, y_max = plot_utils.get_ranges(xs[i_indices, j_indices], ys[i_indices, j_indices])
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)

    #draw a label
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    dx = xmax - xmin
    dy = ymax - ymin
    plt.annotate(label, xy = (xmax - 0.1 * dx, ymax - 0.1 * dy))


    #plot colorbar
    if show_colorbar:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(plot_axes)
        cax = divider.append_axes("right", "8%", pad="3%")

        cb = plt.colorbar(image, ticks = colorbar_tick_locator,
                          orientation = colorbar_orientation,
                          format = colorbar_label_format, drawedges = True,
                          cax=cax
                          )

        cb.outline.set_visible(False)
        cb.ax.set_title(units)

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

    save_to_txt = False

    folder_path = 'data/streamflows/hydrosheds_euler9/'
    i_indices, j_indices = data_select.get_indices_from_file(folder_path + 'aex_discharge_1970_01_01_00_00.nc')
    significance_counter = None
    plt.subplots_adjust(left = 0., hspace = 0.2, wspace = 0.2)

    labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

    ##for querying high flow data for saving to text file
    current_query = None
    future_query = None
    current_highs = None
    future_highs = None
    if level_type == 'high' and save_to_txt:
        high_period_start_month = 3
        high_period_end_month = 7

        current_start_date = datetime(1970,1,1,0,0)
        current_end_date = datetime(1999,12,31,0,0)

        future_start_date = datetime(2041,1,1,0,0)
        future_end_date = datetime(2070,12,31,0,0)

        future_query = QueryObject()
        future_query.start_date = future_start_date
        future_query.end_date = future_end_date
        future_query.event_duration = timedelta(days = 1)
        future_query.start_month = high_period_start_month
        future_query.end_month = high_period_end_month

        current_query = QueryObject()
        current_query.start_date = current_start_date
        current_query.end_date = current_end_date
        current_query.event_duration = timedelta(days = 1)
        current_query.start_month = high_period_start_month
        current_query.end_month = high_period_end_month


    current_id_to_changes = {}
    all_current = []
    all_future = []
    all_stds_current = []
    all_stds_future = []
    for k, current_id in enumerate(members.current_ids):
        if level_type == 'high' and save_to_txt:
            current_path = folder_path + '{0}_discharge_1970_01_01_00_00.nc'.format(current_id)
            future_path = folder_path + '{0}_discharge_2041_01_01_00_00.nc'.format(members.current2future[current_id])
            current_data, times_current, x_indices, y_indices = data_select.get_data_from_file(current_path)
            future_data, times_future, x_indices, y_indices = data_select.get_data_from_file(future_path)

            current_highs = data_select.get_period_maxima_query(current_data, times_current, current_query)
            future_highs = data_select.get_period_maxima_query(future_data, times_future, future_query)

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
        if significance_counter is None:
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

        all_current.append(return_levels_current)
        all_future.append(return_levels_future)

        all_stds_current.append(stdevs_current)
        all_stds_future.append(stdevs_future)

        change[sign_index] /= return_levels_current[sign_index]
        change[sign_index] *= 100.0

        
        plt.subplot(3, 2, k + 1)
        if not (level_type == "high"):
            delta = 100
        else:
            delta = 50


        not_significant = np.ones(change.shape) * 0.75
        not_significant = np.ma.masked_where(condition, not_significant)


        if level_type == 'high':
            assert np.all(return_levels_current > 0)
            assert np.all(stdevs_current >= 0)
            assert np.all(stdevs_future >= 0)

        #temp change to condition
        #change = np.ma.masked_where(np.logical_not(condition), change)
        print 'Plotting: current %s, future %s' % (current_id, future_id)

        current_id_to_changes[current_id] = change

        plot(change , i_indices, j_indices, xs, ys,
                    title = "", label = labels[k],
                    color_map = mycolors.get_red_blue_colormap(ncolors = 20), units = '%',
                    basemap = basemap, minmax = (-delta, delta),
                    colorbar_label_format = '%d',
                    upper_limited = True, colorbar_tick_locator = LinearLocator(numticks = 11),
                    not_significant_mask = None #not_significant
                    )

        if return_period == 10 and level_type == 'high' and save_to_txt:
            txt_saver.save_to_file_rls_and_sign(current_id, return_period,
                                      return_levels_current, return_levels_future,
                                      stdevs_current, stdevs_future,
                                      condition, current_highs, future_highs)
    



        
    plt.subplot(3,2,6)

    plot_sign_count = False
    if plot_sign_count: #if plotting significance count
        significance_counter = np.ma.masked_where(significance_counter == 0, significance_counter)
        plot(significance_counter, i_indices, j_indices, xs, ys,
             title = 'Significance Count', label = labels[5], minmax = (1,6),
             color_map = mycolors.get_sign_count_cmap(ncolors = 5), basemap = basemap,
             colorbar_tick_locator = MaxNLocator(nbins = 5),
             colorbar_label_format = '%d'
             )

        #TODO plot +/-
        plus_change = None
        minus_change = None

        for current_id, the_change in current_id_to_changes.iteritems():
            if plus_change is None:
                plus_change = (the_change > 0)
                minus_change = (the_change < 0)
            else:
                plus_change = np.logical_and(the_change > 0, plus_change)
                minus_change = np.logical_and(the_change < 0, minus_change)

        #should be at least one member with significant changes
        plus_change = np.logical_and(plus_change, significance_counter > 0)
        minus_change = np.logical_and(minus_change, significance_counter > 0)

        x_interest = xs[i_indices, j_indices]
        y_interest = ys[i_indices, j_indices]


        x_plus = x_interest[plus_change]
        y_plus = y_interest[plus_change]
        x_minus = x_interest[minus_change]
        y_minus = y_interest[minus_change]

        basemap.scatter(x_plus, y_plus, marker = "+", color = "m", s = 15, zorder = 5, linewidth = 1)

        if len(x_minus) > 0:
            basemap.scatter(x_minus, y_minus, marker = "d", zorder = 6)
    else:
        #plot ensemble mean
        all_current = np.array( all_current )
        all_future = np.array( all_future )
        all_stds_current = np.array( all_stds_current )
        all_stds_future = np.array( all_stds_future )


        mean_current = np.mean(all_current, axis = 0)
        mean_future = np.mean(all_future, axis = 0)
        mean_stds_current = np.mean( all_stds_current, axis = 0 )
        mean_stds_future = np.mean( all_stds_future, axis = 0 )

        if not level_type == "high":
            delta = 100
        else:
            delta = 50

        not_significant = np.absolute(mean_future - mean_current) <= 1.96 * (mean_stds_current + mean_stds_future)
        not_significant = not_significant.astype(int)
        print " sum(not_significant) = ", np.sum(not_significant)
        not_significant = np.ma.masked_where(~(not_significant == 1), not_significant)
        not_significant *= 0.0


        plot((mean_future - mean_current) / mean_current * 100.0, i_indices, j_indices, xs, ys,
                    title = "", label = labels[-1],
                    color_map = mycolors.get_red_blue_colormap(ncolors = 20), units = '%',
                    basemap = basemap, minmax = (-delta, delta),
                    colorbar_label_format = '%d',
                    upper_limited = True, colorbar_tick_locator = LinearLocator(numticks = 11),
                    not_significant_mask = not_significant
                    )



        pass

    plt.savefig('%d_%s_change_rl.png' % (return_period, level_type), bbox_inches='tight')



def get_column(index, lines, sep = ';'):
    """
    get column index from lines, where columns are separated
    with sep
    """
    selector = lambda x: x.split(sep)[index]
    sel1 = map(selector, lines)
    sel2 = map(float, sel1)
    return sel2


def plot_naveed_troubled_points(path = 'data/data_txt_naveed/TroubledGridCells_AtLeast3Sims.csv'):
    f = open(path)

    lines = f.readlines()
    if len(lines) == 1:
        lines = lines[0].split('\r')
    f.close()

    #skip header
    while len(lines) > 0:
        line = lines.pop(0)
        if 'stationsri' in line.lower():
            break

    print lines

    folder_path = 'data/streamflows/hydrosheds_euler9/'
    i_indices, j_indices = data_select.get_indices_from_file(folder_path + 'aex_discharge_1970_01_01_00_00.nc')


    cols = [1,2,3]
    names = ['C', 'F', 'Both']


    color_map = ListedColormap(['r','b', 'g'])
    plt.subplots_adjust(hspace = 0.0001)
    for c, name in zip(cols, names):
        plt.subplot(2,2,c)
        to_plot = np.ma.masked_all(xs.shape)
        data = get_column(c, lines)
        print data
        for i, j, v in zip(i_indices, j_indices, data):
            if v: #no color for 0
                to_plot[i, j] = v
            else:
                to_plot[i, j] = 3

        basemap.pcolormesh(xs, ys, to_plot.copy(), cmap = color_map,
                          shading = 'flat', rasterized = False )


        plot_basin_boundaries_from_shape(basemap, plotter = plt, linewidth = 1, edge_color = 'k')
        basemap.drawcoastlines(linewidth = 0.5)

        plot_utils.draw_meridians_and_parallels(basemap, step_degrees = 30)

        plt.title(name)

        ymin, ymax = plt.ylim()
        plt.ylim(ymin + 0.05 * (ymax - ymin) , ymax * 0.25)

        xmin, xmax = plt.xlim()
        plt.xlim(xmin + (xmax - xmin) * 0.55, 0.72*xmax)

    cb = plt.colorbar(shrink = 0.5, format = '%d')
    cb.set_ticks([1.25, 2, 2.75])
    cb.set_ticklabels([1,2,3])
    

    plt.savefig('TroubledGridCells_AtLeast3Sims.pdf', bbox_inches = 'tight')

    pass

def plot_naveed_data(path = 'data/data_txt_naveed/3_all_for_plotting.csv'):
    """
    plotting data calculated by Naveed
    """
    f = open(path)

    lines = f.readlines()
    f.close()

    #skip header
    while len(lines) > 0:
        line = lines.pop(0)
        if line.startswith('stationSrI'):
            break

    

    folder_path = 'data/streamflows/hydrosheds_euler9/'
    i_indices, j_indices = data_select.get_indices_from_file(folder_path + 'aex_discharge_1970_01_01_00_00.nc')


    return_periods = [2,10,30]
    level_cols = [4,5,6]
    sig95_cols = [9,10,11]
    sig90_cols = [14, 15, 16]

    sets = zip(level_cols, return_periods, sig95_cols, sig90_cols)

    subplot_count = 1
    for lev_col, period, sig95_col, sig90_col in sets:
        for sig_col in [sig95_col, sig90_col]:
            ret_lev = np.array(get_column(lev_col, lines))
            plt.subplot(3,2,subplot_count)
            significance = 1 - np.array(get_column(sig_col, lines))
            ret_lev = np.ma.masked_where(significance == 1, ret_lev)

            significance *= 0.75
            significance = np.ma.masked_where(significance == 0, significance)

            delta = 50
            sig_level = '95%' if sig_col == sig95_col else '90%'
            plot(ret_lev , i_indices, j_indices, xs, ys,
                        title = 'T = %d year, conf. (%s)' % (period, sig_level),
                        color_map = mycolors.get_red_blue_colormap(ncolors = 16), units = '%',
                        basemap = basemap, minmax = (-delta, delta),
                        colorbar_label_format = '%d',
                        upper_limited = True, colorbar_tick_locator = LinearLocator(numticks = 9),
                        not_significant_mask = significance, show_colorbar = (subplot_count == 6)
                        )
            subplot_count += 1
    plt.savefig('3_all_for_plotting.pdf', bbox_inches = 'tight')
    
    return

    plt.figure()
    b = Basemap(resolution = 'i')
    plt.subplot(2,1,1)
    medians = get_column(1, lines)
    to_plot = np.ma.masked_all(polar_stereographic.lons.shape)
    for i, j, med in zip(i_indices, j_indices, medians):
        to_plot[i,j] = med

    b.pcolormesh(polar_stereographic.lons, polar_stereographic.lats, to_plot)
    b.drawcoastlines()
    plt.colorbar(shrink = 0.5)

    cond = ~to_plot.mask
    min_lon = polar_stereographic.lons[cond].min()
    max_lon = polar_stereographic.lons[cond].max()

    min_lat = polar_stereographic.lats[cond].min()
    max_lat = polar_stereographic.lats[cond].max()

    marginx = 0.05 * (max_lon - min_lon)
    marginy = 0.05 * (max_lat - min_lat)

    plt.xlim(min_lon - marginx, max_lon + marginx)
    plt.ylim(min_lat - marginy, max_lat + marginy)
    plt.title('median change (%)')


    plt.subplot(2,1,2)
    pvalues = get_column(2, lines)
    to_plot = np.ma.masked_all(polar_stereographic.lons.shape)
    for i, j, pvalue in zip(i_indices, j_indices, pvalues):
        to_plot[i,j] = pvalue


    b.pcolormesh(polar_stereographic.lons, polar_stereographic.lats, to_plot)
    b.drawcoastlines()
    plt.colorbar(shrink = 0.5)
    plt.title('p-value for median change')
    marginx = 0.05 * (max_lon - min_lon)
    marginy = 0.05 * (max_lat - min_lat)
    plt.xlim(min_lon - marginx, max_lon + marginx)
    plt.ylim(min_lat - marginy, max_lat + marginy)
    plt.savefig('median.pdf', bbox_inches = 'tight')


    pass

def main():

    high_ret_periods = [10, 30]
    low_ret_periods = [2, 5]


    for ret_period in high_ret_periods:
        calculate_and_plot(ret_period, gevfit.get_high_ret_level_stationary)

    for ret_period in low_ret_periods:
        calculate_and_plot(ret_period, gevfit.get_low_ret_level_stationary)
       

    pass


if __name__ == "__main__":
    #tmp
    #plot_naveed_data()
    #plot_naveed_troubled_points()
    main()
    print "Hello World"
