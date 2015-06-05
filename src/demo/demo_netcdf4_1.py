__author__="huziy"
__date__ ="$9-Mar-2011 10:02:46 AM$"

import demo_netcdf4
import demo_netcdf4_2
import netCDF4 as nc
import application_properties
if __name__ == "__main__":
    print "here"
    application_properties.set_current_directory()
    demo_netcdf4.test()
    demo_netcdf4_2.print_version()
    path = 'data/streamflows/hydrosheds_euler10_spinup100yrs/aex_discharge_1970_01_01_00_00.nc'
    ds = nc.Dataset(path)
    x = ds.variables['water_discharge'][:]
    print x.shape
    ds.close()
    print "Hello World"
