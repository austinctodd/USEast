#===============================================================================
#
# SCRIPT :  plot_transit_times.py
#
# PURPOSE : Ingest input data from US East water age model and plot the variables
#           on a map.  Try masking out non-shelf values
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 16 JANUARY 2015.
#
#===============================================================================

transit_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
               'avg_transit_times.nc'
Artime_file  = '/Volumes/Black_box/Data/USeast-rtime/output/ATL/'+\
               'rtime_ATL.nc'
ATL_file     = '/Volumes/Black_box/Data/USeast/Data/grd/grid_ATLscope.nc'
GOM_file     = '/Volumes/Black_box/Data/USeast/Data/grd/grid_GOM_shelf_scope.nc'
plot_dir     = '/Volumes/Black_box/Data/PLOTS/USeast-age/rtime/All/'

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script plot_transit_times.py"
print "----------------------------------------------------------------------"
print ""

##########################################################################
# USER: DO NOT MAKE CHANGES BELOW (if you do, you're on your own!)
##########################################################################

#-------------------------------------------------------------------------------
# Define all required libraries, routines, and modules
#-------------------------------------------------------------------------------
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import sys
import mpl_toolkits.basemap
import os
#import PIL

#-------------------------------------------------------------------------------
# Define all required sub-libraries, sub-routines, and sub-modules
#-------------------------------------------------------------------------------
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
#from numpy import append, arange, dtype, linspace, meshgrid, ndarray

#-------------------------------------------------------------------------------
# Open grid file and read in variables
#-------------------------------------------------------------------------------
print 'Ingesting data from file %s ...' % (ATL_file)
grid_data=Dataset(ATL_file,mode='r')
hmask  = grid_data.variables['mask_rho' ][:]
Ascope = grid_data.variables['scope_rho'][:]
mask   = grid_data.variables['mask_rho' ][:,:]
lon    = grid_data.variables['lon_rho'  ][:,:]
lat    = grid_data.variables['lat_rho'  ][:,:]
h      = grid_data.variables['h'        ][:,:]
pm     = grid_data.variables['pm'       ][:]
pn     = grid_data.variables['pn'       ][:]
grid_data.close()

scope=Ascope
lon2=lon-0.5*(lon[1,2]-lon[1,1])
lat2=lat-0.5*(lat[2,1]-lat[1,1])

#-------------------------------------------------------------------------------
# Open transit time and fwd file
#-------------------------------------------------------------------------------
transit_data=Dataset(transit_file,mode='r')
Artime_data =Dataset(Artime_file,mode='r')

#-------------------------------------------------------------------------------
# Mask the age variables to the shelf
#-------------------------------------------------------------------------------
scope2=scope+hmask
scope=1-scope
hmask=1-hmask
h=np.ma.array(h,mask=hmask)
cm.gist_ncar.set_bad('white',alpha=0)

#-------------------------------------------------------------------------------
# Open age file
#-------------------------------------------------------------------------------
day_count=[0,31,28,31,30,31,30,31,31,30,31,30,31]
fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor='w')
cm.spectral_r.set_bad('grey',alpha=0)

for i in range(0,1):

  #-----------------------------------------------------------------------------
  # Load and average data for each season
  #-----------------------------------------------------------------------------
  print 'Plotting data for Month %2d' % (i+1)
  stind=sum(day_count[:i+1])*8
  enind=sum(day_count[:i+2])*8
  Ttime=np.zeros(shape=(482,402))
  age   =           transit_data.variables['mean_age_03'][i,35,:,:]
  Artime=np.nanmean(Artime_data.variables['age_02'      ][stind:enind,35,:,:],axis=0)
  
  Artime=Artime*Ascope
  Ttime =age+Artime

  #-----------------------------------------------------------------------------
  # Mask the age variables to the shelf
  #-----------------------------------------------------------------------------
  Ttime=np.ma.array(Ttime,mask=scope)
  #rtime=np.ma.array(rtime,mask=np.isnan(rtime))
  
  #-----------------------------------------------------------------------------
  # Plot Data using Basemap mercator projection
  #-----------------------------------------------------------------------------
  print 'Plotting data.'
  plt.clf()
  m=Basemap(projection='merc',llcrnrlat=26,urcrnrlat=46,\
          llcrnrlon=-81.75,urcrnrlon=lon[481,376],lat_ts=20,resolution='i')
    #m=Basemap(projection='merc',llcrnrlat=18,urcrnrlat=31,\
    #        llcrnrlon=-98.5,urcrnrlon=-80,lat_ts=20,resolution='f')

  # Draw a think border around the whole map
  #m.drawmapboundary(linewidth=3)

  # Plot Data
  xi, yi = m(lon2,lat2)
  cs = m.pcolormesh(xi,yi,Ttime,cmap=cm.spectral_r,vmin=0, vmax = 365.25*5)
  xi, yi = m(lon,lat)
  cM = m.contour(   xi,yi,scope2,[1,1],colors='k',linewidth=2.0)
  cn = m.contour(   xi,yi,Ttime,[182.625,365.25,730.5,913.125,1095.75],colors='0.5',linewidth=1.0)

  # Add Grid Lines
  #m.drawparallels(np.arange(0.,90.,10.), labels=[1,0,0,0],\
  #                fontsize=10, linewidth=0.75, color='.5')
  #m.drawmeridians(np.arange(-105.,-35.,10.), labels=[0,0,0,1],\
  #                fontsize=10, linewidth=0.75, color='.5')

  # Add Coastlines, States, and Country Boundaries
  m.drawcoastlines()
  #m.drawstates()
  #m.drawcountries()
  #m.fillcontinents(color='grey')
  
  # Add Colorbar
#  clabels=np.array([0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3])
  clabels=np.array([0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5])
  clabels=clabels*365.25
  cbar = m.colorbar(cs, location='bottom', pad="4%",ticks=clabels)#,extend='max')
#  cbar.ax.set_xticklabels([0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3])
  cbar.ax.set_xticklabels([0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5])
  cbar.set_label('Transit Time (years)')
  
  # Add Title
  #  plt.title('Season')

  # Save figure to file
  plot_file=plot_dir+'month_'+str(i+1)+'.png'
  print 'Saving figure to file %s ' % plot_file
  plt.savefig('/Users/actodd/Desktop/temp.png',dpi=150,bbox_inches='tight')

transit_data.close()
Artime_data.close()