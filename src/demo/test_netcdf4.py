__author__="huziy"
__date__ ="$Mar 5, 2011 9:09:42 PM$"

#from netCDF4 import Dataset
import test_netcdf4_2

import netCDF4 as nc
import sys

def test():
    path = '/home/huziy/skynet3_rech1/Netbeans Projects/Python/GevFit/data/streamflows/hydrosheds_euler5/aet_discharge_1970_01_01_00_00.nc'
    fpin = nc.Dataset(path).variables['time'][:]
    
    sys.stdout.write('netcdf4-python version: %s\n'%nc.__version__)
    sys.stdout.write('HDF5 lib version:       %s\n'%nc.__hdf5libversion__)
    sys.stdout.write('netcdf lib version:     %s\n'%nc.__netcdf4libversion__)
    pass

if __name__ == "__main__":
    test()
    print "Hello World"
