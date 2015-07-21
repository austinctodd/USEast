#===============================================================================
#
# SCRIPT :  plot_ATL_rtime_transparent.py
#
# PURPOSE : Ingest input data from US East water age model and plot the residence
#           time on a map.  If month 1, plot with US East Coast basemap.
#           Otherwise, plot with transparent background and no coastline.
#
# METHOD :  USER MUST SPECIFY THE DEPTH
#
# HISTORY : Created by Austin Todd on 30 June 2015 for personal use.
#
#===============================================================================

input_file = '/Volumes/Black_box/Data/USeast-rtime/output/ATL/'+\
             'rtime_ATL.nc'
grid_file  = '/Volumes/Black_box/Data/USeast/Data/grd/grid_ATLscope.nc'
plot_dir   = '/Volumes/Black_box/Data/PLOTS/USeast-age/rtime/ATL/bottom/stacked/'

s_layer = 0 # 0=bottom, 35=surface

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
if s_layer ==0:
  sslayer='bottom'
elif s_layer==35:
  sslayer='surface'
else:
  sslayer='RANDOM LAYER'

print "----------------------------------------------------------------------"
print " RUNNING: Script plot_ATL_rtime_transparent.py"
print ""
print " DEPTH : "+sslayer
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
print 'Ingesting data from file %s ...' % (grid_file)
grid_data=Dataset(grid_file,mode='r')
hmask = grid_data.variables['mask_rho' ][:,:401-24]
scope = grid_data.variables['scope_rho'][:,:401-24]
h     = grid_data.variables['h'        ][:,:401-24]
lon   = grid_data.variables['lon_rho'  ][:,:401-24]
lat   = grid_data.variables['lat_rho'  ][:,:401-24]
grid_data.close()

lon2=lon-0.5*(lon[1,2]-lon[1,1])
lat2=lat-0.5*(lat[2,1]-lat[1,1])

#-------------------------------------------------------------------------------
# Open input data file
#-------------------------------------------------------------------------------
input_data=Dataset(input_file,mode='r')

#-------------------------------------------------------------------------------
# Mask the age variables to the shelf
#-------------------------------------------------------------------------------
scope2=scope+hmask
scope=1-scope
hmask=1-hmask
h=np.ma.array(h,mask=hmask)
cm.gist_ncar.set_bad('white',alpha=0.0)

day_count=[0,31,28,31,30,31,30,31,31,30,31,30,31]

for i in range(0,12):

  #-----------------------------------------------------------------------------
  # Load and average data for each season
  #-----------------------------------------------------------------------------
  print 'Loading data for season %1d' % (i+1)
  stind=sum(day_count[:i+1])*8
  enind=sum(day_count[:i+2])*8
  rtime=np.nanmean(input_data.variables['age_02'][stind:enind,s_layer,:,:401-24],axis=0)

  #-----------------------------------------------------------------------------
  # Mask the age variables to the shelf
  #-----------------------------------------------------------------------------
  rtime=np.ma.array(rtime,mask=scope)
  rtime=np.ma.array(rtime,mask=np.isnan(rtime))
  
  # Plot Data
  if i==0:
    #-----------------------------------------------------------------------------
    # Plot Data using Basemap mercator projection
    #-----------------------------------------------------------------------------
    fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor='w')
    plt.clf()
    m=Basemap(projection='merc',llcrnrlat=20,urcrnrlat=lat[481,401-25]+8,\
              llcrnrlon=-90.75,urcrnrlon=lon[481,401-25]+5,lat_ts=20,resolution='f')

    # Draw a thick border around the whole map
    m.drawmapboundary(linewidth=3)
  
    xi, yi = m(lon2,lat2)
    cs = m.pcolormesh(xi,yi,rtime,cmap=cm.gist_ncar,vmin=0, vmax = 365.25*4)
    xi, yi = m(lon,lat)
    cM = m.contour(   xi,yi,scope2,[1,1],colors='k',linewidth=2.0)
    cn = m.contour(   xi,yi,rtime,[365.25,547.875,730.5],colors='0.5',linewidth=1.0)

    # Add Grid Lines
    #m.drawparallels(np.arange(0.,90.,5.), labels=[1,0,0,0],\
    #                fontsize=10, linewidth=0.75, color='.5')
    #m.drawmeridians(np.arange(-105.,-35.,5.), labels=[0,0,0,1],\
    #                fontsize=10, linewidth=0.75, color='.5')

    # Add Coastlines, States, and Country Boundaries
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.fillcontinents(color='grey')

    # Add Colorbar
    clabels=np.array([0,0.5,1,1.5,2,2.5,3,3.5,4])
    clabels=clabels*365.25
    cbar = m.colorbar(cs, location='right', pad="2%",ticks=clabels,extend='max')
    cbar.set_label('Residence Time (years)')
    cbar.ax.set_yticklabels([0,0.5,1,1.5,2,2.5,3,3.5,4])
  else:
    fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor=None)
    plt.clf()
    fig.patch.set_alpha(0.0)
    fig.figurePatch.set_alpha(0.0)
    ax = fig.add_subplot(111)
    ax.patch.set_alpha(0.0)
    m=Basemap(projection='merc',llcrnrlat=26.5,urcrnrlat=lat[481,401-25],\
              llcrnrlon=-81.75,urcrnrlon=lon[481,401-25]+1,lat_ts=20,resolution='i')
              #m.set_alpha(0.0)
    xi, yi = m(lon2,lat2)
    cs = m.pcolormesh(xi,yi,rtime,cmap=cm.gist_ncar,vmin=0, vmax = 365.25*4)
    xi, yi = m(lon,lat)
    m.contour(   xi,yi,scope2,[1,1],colors='k',linewidth=2.0)
    m.contour(   xi,yi,rtime,[365.25,547.875,730.5],colors='0.5',linewidth=1.0)
    cm.gist_ncar.set_bad('white',alpha=0.0)

  # Save figure to file
  plot_file=plot_dir+'test_'+str(i)+'.png'
  print 'Saving figure to file %s ' % plot_file
  plt.savefig(plot_file,dpi=300,bbox_inches='tight',Transparent=True)

