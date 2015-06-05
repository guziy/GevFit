__author__="huziy"
__date__ ="$2 fevr. 2011 16:12:02$"


import matplotlib.pyplot as plt
import numpy as np
import pickle

import members
import gevfit

def calculate_forcing_differences(return_period = 10,
                                  prefix = 'gev_params_stationary',
                                  postfix = ''
                                  ):

    level_fields = []


    file = prefix + '_' + members.control_id + postfix
    pars_set = pickle.load(open(file))

    control_field = np.zeros((len(pars_set)))
    for pos, pars in enumerate(pars_set):
        if 'high' in postfix:
            control_field[pos] = gevfit.get_high_ret_level_stationary(pars, return_period)
        else:
            control_field[pos] = gevfit.get_low_ret_level_stationary(pars, return_period)


    for id in members.current_ids:
        file = prefix + '_' + id + postfix
        pars_set = pickle.load(open(file))

        level_field = np.zeros((len(pars_set)))
        for pos, pars in enumerate(pars_set):
            if 'high' in postfix:
                level_field[pos] = gevfit.get_high_ret_level_stationary(pars, return_period)
            else:
                level_field[pos] = gevfit.get_low_ret_level_stationary(pars, return_period)

        level_fields.append(level_field)


    diff = np.zeros((len(pars_set)))
    for level_field in level_fields:
        diff += (level_field - control_field) # / control_field * 100

    print(np.min(diff), np.max(diff))

    print('min(control_field): ', np.min(control_field))
    return diff / float( len(level_fields) )


def main():
    #calculate and plot boundary forcing errors
    param_file = 'gev_params_stationary'
    years = [10, 30]
    postfixes = ['low', 'high']

    limits = { 10 : 135 , 30 : 770 }

    plt.subplots_adjust(wspace=0.05)
    for i, return_period in enumerate(years):
        for j, postfix in enumerate(postfixes):
            bfe = calculate_forcing_differences(return_period, prefix = param_file, postfix = '_' + postfix)
            plt.subplot( len(postfixes), len(years), i * len(postfixes) + j + 1 )
            gevfit.plot_data(data = bfe, imagefile = None, units = '${\\rm m^3/s}$',
                            minmax = (-limits[return_period], limits[return_period]))
            plt.title( '%s flow, return period %d years' % (postfix, return_period), {'fontsize' : 25} )

    plt.savefig('forcing_uncert.png')


if __name__ == "__main__":
    main()
    print("Hello World")
