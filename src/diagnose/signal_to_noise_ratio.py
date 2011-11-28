from matplotlib.font_manager import FontProperties
import application_properties

__author__="huziy"
__date__ ="$Aug 8, 2011 10:39:54 PM$"


from matplotlib.ticker import LinearLocator
import calculate_significance_from_bootstrap as csfb
import gevfit
import members
import data_select
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

import pylab
import shape.basin_boundaries as basin_boundaries
def apply_plot_params():
    inches_per_pt = 1.0 / 72.27               # Convert pt to inch
    golden_mean = (np.sqrt(5.0) - 1.0) / 2.0       # Aesthetic ratio
    fig_width = 800 * inches_per_pt          # width in inches
    fig_height = fig_width * golden_mean      # height in inches
    fig_size = [fig_width,  1.5 * fig_height]

    font_size = 16

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



def plot_cv_for_return_levels_current_and_future_and_change():
    low_periods = [2, 5]
    high_periods = [10, 30]


    extreme_types = ['high', 'low']
    extreme_type_to_periods = dict(zip(extreme_types, [high_periods, low_periods]))


    i_indices, j_indices = data_select.get_indices_from_file()
    current_cvs = {}
    future_cvs = {}
    for extreme_type in extreme_types:
        current_cvs[extreme_type] = {} #to save return level fields for a given extreme type
        future_cvs[extreme_type] = {}

        for return_period in extreme_type_to_periods[extreme_type]:
            c_rl_fields = np.zeros((len(members.current_ids), len(i_indices)))
            f_rl_fields = np.zeros((len(members.current_ids), len(i_indices)))
            for member_num, the_member in enumerate(members.current_ids):
                pars_current = csfb.get_pars_for_member_and_type(the_member, level_type = extreme_type)
                pars_future = csfb.get_pars_for_member_and_type(members.current2future[the_member], level_type = extreme_type)

                #calculate changes for each pair of members and each position
                npos = len(pars_current)
                for pos, pc, pf in zip(xrange(npos), pars_current, pars_future):
                    if extreme_type == 'high':
                        c = gevfit.get_high_ret_level_stationary(pc, return_period)
                        f = gevfit.get_high_ret_level_stationary(pf, return_period)
                    else:
                        c = gevfit.get_low_ret_level_stationary(pc, return_period)
                        f = gevfit.get_low_ret_level_stationary(pf, return_period)

                    c_rl_fields[member_num, pos] = c
                    f_rl_fields[member_num, pos] = f

                    pass


            current_cvs[extreme_type][return_period] = np.std( c_rl_fields, axis = 0 ) / np.mean( c_rl_fields, axis = 0 )
            future_cvs[extreme_type][return_period] = np.std( f_rl_fields, axis = 0 ) / np.mean( f_rl_fields, axis = 0 )


    #do plotting
    apply_plot_params()
    _plot_map_as_subplots(i_indices, j_indices, current_cvs, title = " Current Climate, CV ")
    plt.savefig('cv_for_return_levels_current.png')
    _plot_map_as_subplots(i_indices, j_indices, future_cvs, title = "Future Climate, CV")
    plt.savefig('cv_for_return_levels_future.png')
    pass






def _plot_map_as_subplots(i_indices, j_indices, extremetype_retperiod_cv_map, title = ""):
    plt.figure()
    i_subplot = 1
    plt.subplots_adjust(hspace = 0.4, wspace = 0.000)
    #TODO: plot CV for the current and future climates

    plt.figtext(0.5, 0.05, title, horizontalalignment='center')

    max_value = 0.1
    for extreme_type, ret_period_to_cv in extremetype_retperiod_cv_map.iteritems():
        #max_value = 0.5 if extreme_type == 'low' else 2.5
        for ret_period, cv_field in ret_period_to_cv.iteritems():
            to_plot = np.ma.masked_all(csfb.xs.shape)

            for i, j, cv in zip(i_indices, j_indices, cv_field):
                to_plot[i, j] = cv


            plt.subplot(2,2, i_subplot)
            i_subplot += 1
            cMap = mpl.cm.get_cmap('jet', 7)
            cMap.set_over(color = '0.5')

            basin_boundaries.plot_basin_boundaries_from_shape(csfb.basemap, plotter = plt, linewidth = 1.6)
            csfb.basemap.pcolormesh(csfb.xs, csfb.ys, to_plot, vmin = 0, vmax = max_value, cmap = cMap)
            # extend choices are "both", "min", "max" and "neither"
            plt.colorbar(extend = 'max', ticks = LinearLocator(numticks = 8), format = "%.2f")
            csfb.basemap.drawcoastlines(linewidth = 0.2)

            plt.title('Return period: {0} years,\n {1} flow event.'.format(ret_period, extreme_type))
            ymin, ymax = plt.ylim()
            plt.ylim(ymin + 0.05 * (ymax - ymin) , ymax * 0.25)

            xmin, xmax = plt.xlim()
            plt.xlim(xmin + (xmax - xmin) * 0.55, 0.72*xmax)


def plot_signal_to_noise_ratio():
    low_periods = [2, 5]
    high_periods = [10, 30]
    apply_plot_params()

    extreme_types = ['high', 'low']
    extreme_type_to_periods = dict(zip(extreme_types, [high_periods, low_periods]))

    i_indices, j_indices = data_select.get_indices_from_file()
    i_subplot = 1

    plt.subplots_adjust(hspace = 0.4, wspace = 0.000)
    for extreme_type in extreme_types:
        for return_period in extreme_type_to_periods[extreme_type]:

            changes = np.zeros((len(members.current_ids), len(i_indices)))
            for member_num, the_member in enumerate(members.current_ids):

                pars_current = csfb.get_pars_for_member_and_type(the_member, level_type = extreme_type)
                pars_future = csfb.get_pars_for_member_and_type(members.current2future[the_member], level_type = extreme_type)

                
                #calculate changes for each pair of members and each position
                npos = len(pars_current)
                for pos, pc, pf in zip(xrange(npos), pars_current, pars_future):
                    if extreme_type == 'high':
                        c = gevfit.get_high_ret_level_stationary(pc, return_period)
                        f = gevfit.get_high_ret_level_stationary(pf, return_period)
                    else:
                        c = gevfit.get_low_ret_level_stationary(pc, return_period)
                        f = gevfit.get_low_ret_level_stationary(pf, return_period)

                    changes[member_num, pos] = f - c
                    pass

            #calculate mean and stdev of the obtained changes
            the_mean = np.mean(changes, axis = 0)
            the_std = np.std(changes, axis = 0)
            print extreme_type
            print 'std = ', the_std.min(), the_std.max()
            print 'mean = ', the_mean.min(), the_mean.max()
            print 20 * '--'

            the_values = np.zeros(the_mean.shape)
            delta = 0

            #change if you want signal to noise ratio, currently it is cv (coefficient of variation 1/(signal-to-noise-ratio))
            the_values[the_std > delta] = the_std[the_std > delta] / np.abs(the_mean[the_std > delta])
            print the_values.min(), the_values.max()
            #the_values = the_mean
            to_plot = np.ma.masked_all(csfb.xs.shape)

            max_value = 0.5
            for i, j, value, dev in zip(i_indices, j_indices, the_values, the_std):
                to_plot[i, j] = value


            #shaded = np.ma.masked_where(to_plot != 0, shaded)
            #to_plot = np.ma.masked_where(to_plot == 0, to_plot)
            print np.any(to_plot == 0)

            plt.subplot(2,2, i_subplot)
            i_subplot += 1
            print 'just before plotting'
            cMap = mpl.cm.get_cmap('jet', 7)
            cMap.set_over(color = '0.5')

            
            basin_boundaries.plot_basin_boundaries_from_shape(csfb.basemap, plotter = plt, linewidth = 1.6)
            csfb.basemap.pcolormesh(csfb.xs, csfb.ys, to_plot, vmin = 0, vmax = max_value, cmap = cMap)
            # extend choices are "both", "min", "max" and "neither"
            plt.colorbar(extend = 'max', ticks = LinearLocator(numticks = 8), format = "%.2f")
            csfb.basemap.drawcoastlines(linewidth = 0.2)
            plt.title('Return period: {0} years,\n {1} flow event.'.format(return_period, extreme_type))
            ymin, ymax = plt.ylim()
            plt.ylim(ymin + 0.05 * (ymax - ymin) , ymax * 0.25)

            xmin, xmax = plt.xlim()
            plt.xlim(xmin + (xmax - xmin) * 0.55, 0.72*xmax)


    #plt.show()

    plt.savefig('cv_for_changes.pdf', bbox_inches='tight')
 


if __name__ == "__main__":
    application_properties.set_current_directory()
    plot_cv_for_return_levels_current_and_future_and_change()
    #plot_signal_to_noise_ratio()
    print "Hello World"
