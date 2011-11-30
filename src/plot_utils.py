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




if __name__ == "__main__":
    print "Hello World"
