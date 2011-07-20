__author__="huziy"
__date__ ="$Apr 1, 2011 7:48:56 PM$"

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

from numpy.random import rand

import mpl_toolkits.basemap as bm


if __name__ == "__main__":

    print matplotlib.__version__
    print matplotlib.__revision__
    print matplotlib.__date__


    print 'basemap'
    print bm.__version__

    print 'numpy'
    print np.__version__


    x = rand(15000,10000)
    #print x
    plt.pcolormesh(x)
    plt.savefig('test_pcolor.png', bbox_inches = 'tight')
    print "Hello World"
