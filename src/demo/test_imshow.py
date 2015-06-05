# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="huziy"
__date__ ="$Mar 2, 2011 9:14:21 AM$"

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def test():
    to_plot = np.ma.masked_all((100, 100))
    for i in range(25,75):
        for j in range(25,75):
            to_plot[i, j] = i ** 2 + j ** 2



    alpha =  np.where(to_plot.mask, 0, 1)
    print(alpha)

    plt.pcolormesh(to_plot, cmap =  mpl.cm.get_cmap('Reds', 10), alpha = None)
    #plt.imshow(to_plot, interpolation = 'hanning', cmap =  mpl.cm.get_cmap('Reds', 10), alpha = alpha)
    plt.show()
    pass


if __name__ == "__main__":
    test()
    print("Hello World")
