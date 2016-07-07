# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 15:38:12 2016

@author: raissaphilibert
"""


import datetime as dt  # Python standard library datetime  module
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from netcdftime import utime
from datetime import  datetime
from repel_labels import *

file = '/Users/raissaphilibert/Dropbox/oceanography_mac/nitrification/write-upsandpresentations/Paper_Benguela_nitrification/WSAfrica_MODIS_sst_20111129.nc4'
'''
NAME
    NetCDF with Python
PURPOSE
    To demonstrate how to read and write data with NetCDF files using
    a NetCDF file from the NCEP/NCAR Reanalysis.
    Plotting using Matplotlib and Basemap is also shown.
PROGRAMMER(S)
    Chris Slocum
REVISION HISTORY
    20140320 -- Initial version created and posted online
    20140722 -- Added basic error handling to ncdump
                Thanks to K.-Michael Aye for highlighting the issue
REFERENCES
    netcdf4-python -- http://code.google.com/p/netcdf4-python/
    NCEP/NCAR Reanalysis -- Kalnay et al. 1996
        http://dx.doi.org/10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2
''''''
NAME
    NetCDF with Python
PURPOSE
    To demonstrate how to read and write data with NetCDF files using
    a NetCDF file from the NCEP/NCAR Reanalysis.
    Plotting using Matplotlib and Basemap is also shown.
PROGRAMMER(S)
    Chris Slocum
REVISION HISTORY
    20140320 -- Initial version created and posted online
    20140722 -- Added basic error handling to ncdump
                Thanks to K.-Michael Aye for highlighting the issue
REFERENCES
    netcdf4-python -- http://code.google.com/p/netcdf4-python/
    NCEP/NCAR Reanalysis -- Kalnay et al. 1996
        http://dx.doi.org/10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2
'''

def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print "\t\ttype:", repr(nc_fid.variables[key].dtype)
            for ncattr in nc_fid.variables[key].ncattrs():
                print '\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr))
        except KeyError:
            print "\t\tWARNING: %s does not contain variable attributes" % key

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print "NetCDF Global Attributes:"
        for nc_attr in nc_attrs:
            print '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print "NetCDF dimension information:"
        for dim in nc_dims:
            print "\tName:", dim 
            print "\t\tsize:", len(nc_fid.dimensions[dim])
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print "NetCDF variable information:"
        for var in nc_vars:
            if var not in nc_dims:
                print '\tName:', var
                print "\t\tdimensions:", nc_fid.variables[var].dimensions
                print "\t\tsize:", nc_fid.variables[var].size
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars
    
nc_f = file  # Your filename
nc_fid = Dataset(nc_f, 'r')  # Dataset is the class behavior to open the file
                             # and create an instance of the ncCDF4 class
nc_attrs, nc_dims, nc_vars = ncdump(nc_fid)
# Extract data from NetCDF file
lats = nc_fid.variables['lat'][:]  # extract/copy the data
lons = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][:]
time_units= nc_fid.variables['time'].units
sst = nc_fid.variables['sst'][:]  # shape is time, lat, lon as shown above

#lon, lat = np.meshgrid(lons, lats)
## HERER WILL TRY CALCULATE A CLIMATOLOGY
#mean_adt=[]
#for i in adt:
#    time_idx = i
#    adt_cyclic, lons_cyclic = addcyclic(adt[time_idx, :, :], lons)
## Create 2D lat/lon arrays for Basemap
#    lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
## Transforms lat/lon into plotting coordinates for projection
#    x, y = m(lon2d, lat2d)
#    gridded[i]=(lon2d,lat2d,adt_cyclic)
#
#    


###gets the unit of time from the netcdf file
cdftime = utime(time_units)
date1 = cdftime.num2date(time)   # turns the times in date (2014,12,6,0,0) 

date = [datet.date() for datet in date1]
date2 = [datet.toordinal() for datet in date]
##date I want YYYY,M,DAY
date_select = ((datetime(2011,11,29)).date()).toordinal()

time_idx =date2.index(date_select)

 # Plot of global temperature on our random day
fig = plt.figure()
fig.subplots_adjust(left=0., right=1., bottom=0., top=0.9)
# Setup the map. See http://matplotlib.org/basemap/users/mapsetup.html
# for other projections.

# Make the plot continuous

m=Basemap(projection='cyl',llcrnrlon=13,llcrnrlat=-36,urcrnrlon=21.,urcrnrlat=-27.,
          #lat_0=-35, lon_0=18,
          resolution='l', area_thresh=1000.0)

m.drawcoastlines()
m.drawmapboundary()
adt_cyclic, lons_cyclic = addcyclic(sst[time_idx, :, :], lons)
# Shift the grid so lons go from -180 to 180 instead of 0 to 360.
#adt_cyclic, lons_cyclic = shiftgrid(180., adt_cyclic, lons_cyclic, start=False)
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
# Transforms lat/lon into plotting coordinates for projection
x, y = m(lon2d, lat2d)
# Plot of air temperature with 11 contour intervals
cs = m.contourf(x, y, adt_cyclic, 11, cmap=plt.cm.viridis)


# plot the mean frontal positions
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
#m.contour(x,y,adt_cyclic,[-1.244],colors='r',linewidth=1.3)# % Sbdy
#m.contour(x,y,adt_cyclic,[-0.943],colors='r',linewidth=1.3)#% SACCF_N

cbar = plt.colorbar(cs, orientation='horizontal', shrink=0.5)
cbar.set_label("%s (%s)" % (nc_fid.variables['sst'].standard_name,\
                          nc_fid.variables['sst'].units))
cur_time=(datetime.fromordinal((date_select)).date()).strftime('%m/%d/%Y')

m.drawcoastlines()
m.drawcountries()
m.fillcontinents(color='coral')
m.drawmapboundary()
parallels = np.arange(-36,-27,5.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[True,False,False,True])
meridians = np.arange(-180,180,5.)
m.drawmeridians(meridians,labels=[True,True,False,True])



plt.title("%s" % (cur_time))
 

#data_lat (samples lats) 
data_lat=[-32]
data_lon=[18]
name =['Station']

x, y = m(data_lon, data_lat)
xy=zip(data_lon,data_lat)
m.scatter(x, y, marker='D',color='m')
#text(data_lon,data_lat,'Station')
                           
fig.savefig('sst_sampple_pos.eps', bbox_inches='tight')
nc_fid.close()