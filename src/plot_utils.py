# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="huziy"
__date__ ="$24 juil. 2010 16:45:38$"

import numpy as np

def draw_meridians_and_parallels(the_basemap, step_degrees = 5.0):
    meridians = np.arange(-180,180, step_degrees)
    parallels = np.arange(-90,90, step_degrees)
    the_basemap.drawparallels(parallels,labels=[0,0,0,0],fontsize=16,linewidth=0.25)
    the_basemap.drawmeridians(meridians,labels=[0,0,0,0],fontsize=16,linewidth=0.25)

def get_ranges(x_interest, y_interest):
    """
    Get region of zoom for a given map
    """
    x_min, x_max = np.min( x_interest ), np.max( x_interest )
    y_min, y_max = np.min( y_interest ), np.max( y_interest )
    dx = 0.1 * ( x_max - x_min )
    dy = 0.1 * ( y_max - y_min )
    return x_min - dx, x_max + dx, y_min - dy, y_max + dy



def get_closest_tick_value(nticks, lower_limit):
    """
    nticks - number of ticks in the colorbar
    lower_limit - is the lower limit of the data to plot [0..1]
    """

    assert 0 <= lower_limit <= 1
    d = 1.0 / float( nticks - 1.0 )
    assert d > 0

    tick_value = 0
    while tick_value <= 1:
        if tick_value <= lower_limit <= tick_value + d:
            if lower_limit - tick_value < tick_value + d - lower_limit:
                return tick_value
            else:
                return tick_value + d
        tick_value += d


def apply_plot_params(font_size = 20, width_pt = 1000, aspect_ratio = 1):
    """
    aspect_ratio = height / (width * golden_mean)
    """
    import pylab
    import math

    if width_pt is not None:
        inches_per_pt = 1.0 / 72.27               # Convert pt to inch
        golden_mean = (math.sqrt(5.0) - 1.0) / 2.0       # Aesthetic ratio
        fig_width = width_pt * inches_per_pt          # width in inches
        fig_height = fig_width * golden_mean      # height in inches
        fig_size = [fig_width,  aspect_ratio * fig_height]
    else:
        inches_per_cm = 1.0 / 2.54
        width_cm = 16.0
        height_cm = 23.0
        fig_size = [ width_cm * inches_per_cm, height_cm * inches_per_cm ]

    params = {
        'axes.labelsize': font_size,
        'font.size':font_size,
        'text.fontsize': font_size,
        'legend.fontsize': font_size,
        'xtick.labelsize': font_size,
        'ytick.labelsize': font_size,
        'figure.figsize': fig_size,
        "axes.titlesize" : font_size
        }

    pylab.rcParams.update(params)


if __name__ == "__main__":
    print "Hello World"
