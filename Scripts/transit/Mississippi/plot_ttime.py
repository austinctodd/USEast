#===============================================================================
#
# SCRIPT :  plot_GOM_shelf_rtime_depav.py
#
# PURPOSE : Ingest input data from US East water age model and plot the depth-
#           averaged residence time on a map of just Gulf of Mexico shelves.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 30 June 2015 for personal use.
#
#===============================================================================

age_file   = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
             'age_pycnocline_avg.nc'
rtime_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
             'rtime_pycnocline.nc'
GOM_file  = '/Volumes/Black_box/Data/USeast/Data/grd/grid_GOM_shelf_scope.nc'
plot_dir   = '/Volumes/Black_box/Data/PLOTS/USeast-age/transit/pcline/Mississippi/'

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script plot_GOM_shelf_rtime_depav.py"
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
print 'Ingesting data from file %s ...' % (GOM_file)
grid_data=Dataset(GOM_file,mode='r')
hmask = grid_data.variables['mask_rho' ][:]
Gscope = grid_data.variables['scope_rho'][:]
h     = grid_data.variables['h'        ][:]
lon   = grid_data.variables['lon_rho'  ][:]
lat   = grid_data.variables['lat_rho'  ][:]
grid_data.close()

scope=Gscope

lon2=lon-0.5*(lon[1,2]-lon[1,1])
lat2=lat-0.5*(lat[2,1]-lat[1,1])

#-------------------------------------------------------------------------------
# Open input data files
#-------------------------------------------------------------------------------
age_data  =Dataset(  age_file,mode='r')
rtime_data=Dataset(rtime_file,mode='r')

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

for i in range(0,12):

  #-----------------------------------------------------------------------------
  # Load and average rtime data for each season
  #-----------------------------------------------------------------------------
  print 'Plotting data for Month %2d' % (i+1)
  stind=sum(day_count[:i+1])*8
  enind=sum(day_count[:i+2])*8
  age  =np.nanmean(  age_data.variables['mean_age_02'][stind:enind,:,:],axis=0)
  rtime=np.nanmean(rtime_data.variables['GOM_shelf_rtime'][stind:enind,:,:],axis=0)
  
  ttime=age+rtime
  
  #-----------------------------------------------------------------------------
  # Mask the age variables to the shelf
  #-----------------------------------------------------------------------------
  ttime=np.ma.array(ttime,mask=scope)
  
  #-----------------------------------------------------------------------------
  # Plot Data using Basemap mercator projection
  #-----------------------------------------------------------------------------
  print 'Plotting data.'
  plt.clf()
  m=Basemap(projection='merc',llcrnrlat=18,urcrnrlat=31,\
            llcrnrlon=-98.5,urcrnrlon=-80,lat_ts=20,resolution='f')

  # Draw a think border around the whole map
  m.drawmapboundary(linewidth=3)

  # Plot Data
  xi, yi = m(lon2,lat2)
  cs = m.pcolormesh(xi,yi,ttime,cmap=cm.gist_ncar,vmin=0, vmax = 365.25*4)
  xi, yi = m(lon,lat)
  cM = m.contour(   xi,yi,scope2,[1,1],colors='k',linewidth=2.0)
  cn = m.contour(   xi,yi,ttime,[182.625,365.25,730.5],colors='0.5',linewidth=1.0)

  # Add Grid Lines
  m.drawparallels(np.arange(0.,90.,10.), labels=[1,0,0,0],\
                  fontsize=10, linewidth=0.75, color='.5')
  m.drawmeridians(np.arange(-105.,-35.,10.), labels=[0,0,0,1],\
                  fontsize=10, linewidth=0.75, color='.5')

  # Add Coastlines, States, and Country Boundaries
  m.drawcoastlines()
  m.drawstates()
  m.drawcountries()
  m.fillcontinents(color='grey')
  
  # Add Colorbar
  clabels=np.array([0,0.5,1,1.5,2,2.5,3,3.5,4])
  clabels=clabels*365.25
  cbar = m.colorbar(cs, location='bottom', pad="4%",ticks=clabels,extend='max')
  cbar.set_label('Transit Time (years)')
  cbar.ax.set_xticklabels([0,0.5,1,1.5,2,2.5,3,3.5,4])

  # Add Title
  #  plt.title('Season')

  # Save figure to file
  plot_file=plot_dir+'month_'+str(i+1)+'.png'
  #plot_file=plot_dir+'annual_mean.png'
  print 'Saving figure to file %s ' % plot_file
  plt.savefig(plot_file,dpi=300,bbox_inches='tight')

age_data.close()
rtime_data.close()