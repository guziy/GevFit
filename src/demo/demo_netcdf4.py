__author__="huziy"
__date__ ="$Mar 5, 2011 9:09:42 PM$"

from netCDF4 import Dataset
import netCDF4 as nc
import sys

import application_properties
import demo_netcdf4_2

def test():

    demo_netcdf4_2.print_version()
    path = 'data/streamflows/output_2d/data1/aex_discharge_1961_01_01_00_00.nc'
    fpin = Dataset(path)
    print fpin.variables.keys()
    x = fpin.variables['time'][:]
    print x.shape
    fpin.close()
    
    sys.stdout.write('netcdf4-python version: %s\n'%nc.__version__)
    sys.stdout.write('HDF5 lib version:       %s\n'%nc.__hdf5libversion__)
    sys.stdout.write('netcdf lib version:     %s\n'%nc.__netcdf4libversion__)
    pass

if __name__ == "__main__":
    application_properties.set_current_directory()
    test()
    print "Hello World"
