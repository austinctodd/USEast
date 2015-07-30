#===============================================================================
#
# SCRIPT :  plot_MS_trajectories.py
#
# PURPOSE : Ingest input data from LTRANS output and plot on a map.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 17 March 2015.
#
#===============================================================================

input_dir  = '/Volumes/Black_box/Data/LTRANS/output/Mississippi/'
GOM_file   = '/Volumes/Black_box/Data/USeast/Data/grd/grid_GOM_shelf_scope.nc'
plot_dir   = '/Volumes/Black_box/Data/PLOTS/LTRANS/Mississippi/frames/'

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script plot_MS_transit_times.py"
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
hmask  = grid_data.variables['mask_rho' ][:]
Gscope = grid_data.variables['scope_rho'][:]
lon    = grid_data.variables['lon_rho'  ][:,:]
lat    = grid_data.variables['lat_rho'  ][:,:]
grid_data.close()

scope=Gscope
lon2=lon-0.5*(lon[1,2]-lon[1,1])
lat2=lat-0.5*(lat[2,1]-lat[1,1])

#-------------------------------------------------------------------------------
# Open and read transit time data
#-------------------------------------------------------------------------------
file1=Dataset(input_dir+'AR_v2.txt',mode='r')
file2=Dataset(input_dir+'MS1_v2.txt',mode='r')
file3=Dataset(input_dir+'MS2_v2.txt',mode='r')
file4=Dataset(input_dir+'MS3_v2.txt',mode='r')
file5=Dataset(input_dir+'MS4_v2.txt',mode='r')
file6=Dataset(input_dir+'MS5_v2.txt',mode='r')

#-------------------------------------------------------------------------------
# Mask the age variables to the shelf
#-------------------------------------------------------------------------------
scope2=scope+hmask
scope=1-scope
hmask=1-hmask
cm.gist_ncar.set_bad('white',alpha=0)

#-------------------------------------------------------------------------------
# Open age file
#-------------------------------------------------------------------------------
fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor='w')
cm.gist_ncar.set_bad('w')

for i in range(0,366):

  #-----------------------------------------------------------------------------
  # Load and average data for each season
  #-----------------------------------------------------------------------------
  #print 'Plotting data for Month %2d' % (i+1)
  Ttime=expos_time*Gscope[99:260,0:200]
  Ttime=Ttime/expos_part

  #-----------------------------------------------------------------------------
  # Mask the age variables to the shelf
  #-----------------------------------------------------------------------------
  Ttime=np.ma.array(Ttime,mask=scope[99:260,0:200])
  Ttime=np.ma.masked_where(expos_part<10,Ttime) #<-- Mask where less than 0.5%
  
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
  xi, yi = m(lon2[99:260,0:200],lat2[99:260,0:200])
  cs = m.pcolormesh(xi,yi,Ttime,cmap=cm.gist_ncar,vmin=0, vmax = 365.25*2.5)
  xi, yi = m(lon,lat)
  cM = m.contour(   xi,yi,scope2,[1,1],colors='k',linewidth=2.0)
  xi, yi = m(lon[99:260,0:200],lat[99:260,0:200])
  cn = m.contour(   xi,yi,Ttime,[182.625,365.25,730.5,913.125,1095.75],colors='0.5',linewidth=1.0)

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
  clabels=np.array([0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5])#,2.75,3])
  clabels=clabels*365.25
  cbar = m.colorbar(cs, location='bottom', pad="4%",ticks=clabels,extend='max')
  cbar.ax.set_xticklabels([0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5])#,2.75,3])
  cbar.set_label('Transit Time (years)')
  
  # Add Title
  #  plt.title('Season')

  # Save figure to file
  #plot_file=plot_dir+'month_'+str(i+1)+'.png'
  plot_file=plot_dir+'annual_mean_expos.png'
  print 'Saving figure to file %s ' % plot_file
  plt.savefig(plot_file,dpi=150,bbox_inches='tight')

