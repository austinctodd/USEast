#===============================================================================
#
# SCRIPT :  plot_monthly_age.py
#
# PURPOSE : Ingest input data from US East water age model and plot the variables
#           on a map.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 2 June 2015.
#
#===============================================================================

input_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
             'age_pycnocline_avg.nc'
GOM_file   = '/Volumes/Black_box/Data/USeast/Data/grd/grid_GOM_shelf_scope.nc'
plot_dir   = '/Volumes/Black_box/Data/PLOTS/USeast-age/age/pcline/GOM/'

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script plot_monthly_age.py"
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
lon   = grid_data.variables['lon_rho'  ][:]
lat   = grid_data.variables['lat_rho'  ][:]
h     = grid_data.variables['h'        ][:]
grid_data.close()

# Make shifted lat/lon values for pcolor to plot center box on lat/lon
lon2=lon-0.5*(lon[1,2]-lon[1,1])
lat2=lat-0.5*(lat[2,1]-lat[1,1])

#-------------------------------------------------------------------------------
# Open transit time and fwd file
#-------------------------------------------------------------------------------
input_data=Dataset(input_file,mode='r')

#-------------------------------------------------------------------------------
# Mask the age variables to the shelf
#-------------------------------------------------------------------------------
hmask=1-hmask
h=np.ma.array(h,mask=hmask)
cm.jet.set_bad('white',alpha=0)

#-------------------------------------------------------------------------------
# Open age file
#-------------------------------------------------------------------------------
day_count=[0,31,28,31,30,31,30,31,31,30,31,30,31]
fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor='w')
cm.gist_ncar.set_bad('grey',alpha=0)

for tdim in range(0,12):

  #------------------------------------------------------------------------
  # Set index limits for current month
  #------------------------------------------------------------------------
  stind = (sum(day_count[0:tdim+1]))*8
  enind = (sum(day_count[0:tdim+2]))*8

  #-----------------------------------------------------------------------------
  # Load and average data for each season
  #-----------------------------------------------------------------------------
  print 'Plotting data for Month %2d' % (tdim+1)
  age   =np.nanmean(input_data.variables['mean_age_03'][stind:enind,:,:],axis=0)
  
  #-----------------------------------------------------------------------------
  # Mask the values over land
  #-----------------------------------------------------------------------------
  age=np.ma.array(age,mask=hmask)
  
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
  cs = m.pcolormesh(xi,yi,age,cmap=cm.gist_ncar,vmin=0, vmax = 365.25*3)
  xi, yi = m(lon,lat)
  #  cM = m.contour(   xi,yi,scope2,[1,1],colors='k',linewidth=2.0)
  cn = m.contour(   xi,yi,age,[182.625,365.25,730.5,913.125],colors='0.5',linewidth=1.0)

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
  clabels=np.array([0,0.5,1,1.5,2,2.5,3])
  clabels=clabels*365.25
  cbar = m.colorbar(cs, location='bottom',pad="4%",ticks=clabels,extend='max')
  cbar.ax.set_xticklabels([0,0.5,1,1.5,2,2.5,3])
  cbar.set_label('Mean age (years)')
  
  # Add Title
  #  plt.title('Season')

  # Save figure to file
  plot_file=plot_dir+'month_'+str(tdim+1)+'.png'
  print 'Saving figure to file %s ' % plot_file
  plt.savefig(plot_file,dpi=150,bbox_inches='tight')

input_data.close()