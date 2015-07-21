#===============================================================================
#
# SCRIPT :  plot_age_panels.py
#
# PURPOSE : Ingest input data from US East water age model and plot the variables
#           on a map.  Try masking out non-shelf values
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 16 JANUARY 2015.
#
#===============================================================================

age_file    = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
              'avg_3hrly.nc'
pcline_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
              'age_pycnocline_avg.nc'
grid_file   = '/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc'
river_file  = '/Volumes/Black_box/Data/USeast/Data/frc/USeast-riverNOVN-12.nc'
plot_dir    = '/Volumes/Black_box/Data/PLOTS/USeast-age/age/pcline/age03/'

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script plot_age_panels.py"
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
import PIL

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
print 'Ingesting data from file %s ...' % (pcline_file)
grid_data=Dataset(grid_file,mode='r')
mask=grid_data.variables['mask_rho' ][:,0:401-25]
lon =grid_data.variables['lon_rho'  ][:,0:401-25]
lonu=grid_data.variables['lon_u'    ][:]
lonv=grid_data.variables['lon_v'    ][:]
lat =grid_data.variables['lat_rho'  ][:,0:401-25]
latu=grid_data.variables['lat_u'    ][:]
latv=grid_data.variables['lat_v'    ][:]
grid_data.close()
mask=1.0-mask

#-------------------------------------------------------------------------------
# Open age file
#-------------------------------------------------------------------------------
pcline_data=Dataset(pcline_file,mode='r')


day_count=[0,31,28,31,30,31,30,31,31,30,31,30,31]
fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor='w')
cm.gist_ncar.set_bad('grey',alpha=0)

for i in range(0,1):

  #-----------------------------------------------------------------------------
  # Load and average data for each season
  #-----------------------------------------------------------------------------
  print 'Loading data for season %1d' % (i+1)
  stind = (sum(day_count[0:(i  )*3+1]))*8
  enind = (sum(day_count[0:(i+1)*3+1]))*8-1
  tmp=pcline_data.variables['mean_age_03'][:,:,0:401-25]
  age=np.nanmean(tmp,axis=0)

  #-----------------------------------------------------------------------------
  # Mask the age variables to the shelf
  #-----------------------------------------------------------------------------
  age=np.ma.array(age,mask=mask)
  age=np.ma.array(age,mask=np.isnan(age))
  
  #-----------------------------------------------------------------------------
  # Plot Data using Basemap mercator projection
  #-----------------------------------------------------------------------------
  print 'Plotting data.'
  plt.clf()
  m=Basemap(projection='merc',llcrnrlat=lat.min(),urcrnrlat=lat.max(),\
            llcrnrlon=lon.min(),urcrnrlon=lon.max(),lat_ts=20,resolution='i')

  # If Lambert Conformal projection preferred (looks similar to Mercator)
  #m=Basemap(width=3850000, height=4600000,
  #          rsphere=(6378137.00,6356752.3142),projection='lcc',resolution='i',\
  #          area_thresh=100., lat_1=5.,lat_2=10., lat_0=lat[230,195],\
  #          lon_0=lon[230,195])

  # Draw a think border around the whole map
  m.drawmapboundary(linewidth=3)

  # Plot Data
  xi, yi = m(lon,lat)
  cs = m.pcolormesh(xi,yi,age,cmap=cm.gist_ncar,vmin=0, vmax = 365.25*2.5)
  cn = m.contour(   xi,yi,age,[365.25,547.875,730.5],colors='0.5',linewidth=1.0)

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
  clabels=np.array([0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5])
  clabels=clabels*365.25
  cbar = m.colorbar(cs, location='right', pad="2%",ticks=clabels,extend='max')
  cbar.set_label('Mean Age (years)')
  cbar.ax.set_yticklabels([0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5])

  # Add Title
  #  plt.title('Season')

  # Save figure to file
  plot_file=plot_dir+'season_'+str(i+1)+'.png'
  plot_file=plot_dir+'annual_mean.png'
  print 'Saving figure to file %s ' % plot_file
  plt.savefig(plot_file,dpi=300,bbox_inches='tight')

