# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="huziy"
__date__ ="$9-Mar-2011 10:02:46 AM$"

import test_netcdf4
import netCDF4 as nc
if __name__ == "__main__":
    print "here"
    test_netcdf4.test()
    path = '/home/huziy/skynet3_rech1/Netbeans Projects/Python/GevFit/data/streamflows/hydrosheds_euler5/aet_discharge_1970_01_01_00_00.nc'
    nc.Dataset(path).variables['water_discharge'][:]
    print "Hello World"
