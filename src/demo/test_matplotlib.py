__author__="huziy"
__date__ ="$Apr 1, 2011 7:48:56 PM$"

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

from numpy.random import rand

import mpl_toolkits.basemap as bm

from map_parameters import polar_stereographic


if __name__ == "__main__":

    print("matplotlib version")
    print(matplotlib.__version__)


    print('basemap')
    print(bm.__version__)

    print('numpy')
    print(np.__version__)


    x = rand(150,100)
    #print x
    plt.pcolormesh(x)
    plt.savefig('test_pcolor.png', bbox_inches = 'tight')
    print("Hello World")


    plt.figure()
    x = polar_stereographic.xs
    y = polar_stereographic.ys
    plt.annotate("+", xy = (x[5, 5],y[5,5]))
    #plt.scatter(x[1:10, 1:10],y[1:10,1:10], marker=r'$_$')
    plt.show()