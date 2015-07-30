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
plot_dir   = '/Volumes/Black_box/Data/PLOTS/LTRANS/Mississippi/'

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
mask   = grid_data.variables['mask_rho' ][:,:]
lon    = grid_data.variables['lon_rho'  ][:,:]
lat    = grid_data.variables['lat_rho'  ][:,:]
grid_data.close()

scope=Gscope
lon2=lon-0.5*(lon[1,2]-lon[1,1])
lat2=lat-0.5*(lat[2,1]-lat[1,1])

#-------------------------------------------------------------------------------
# Open and read transit time data
#-------------------------------------------------------------------------------
file1=input_dir+'AR_v2.txt'
file2=input_dir+'MS1_v2.txt'
file3=input_dir+'MS2_v2.txt'
file4=input_dir+'MS3_v2.txt'
file5=input_dir+'MS4_v2.txt'
file6=input_dir+'MS5_v2.txt'

trans_time=np.zeros(shape=(161,200));
trans_part=np.zeros(shape=(161,200));
expos_time=np.zeros(shape=(161,200));
expos_part=np.zeros(shape=(161,200));

for f in range(0,6):
  if f==0:
    ii,jj,tt,tp,et,ep=np.loadtxt(file1,unpack=True)
  elif f==1:
    ii,jj,tt,tp,et,ep=np.loadtxt(file2,unpack=True)
  elif f==2:
    ii,jj,tt,tp,et,ep=np.loadtxt(file3,unpack=True)
  elif f==3:
    ii,jj,tt,tp,et,ep=np.loadtxt(file4,unpack=True)
  elif f==4:
    ii,jj,tt,tp,et,ep=np.loadtxt(file5,unpack=True)
  elif f==5:
    ii,jj,tt,tp,et,ep=np.loadtxt(file6,unpack=True)

  for i in range(0,ii.size):
    trans_time[jj[i]-100,ii[i]-1]=trans_time[jj[i]-100,ii[i]-1]+tt[i]
    trans_part[jj[i]-100,ii[i]-1]=trans_part[jj[i]-100,ii[i]-1]+tp[i]
    expos_time[jj[i]-100,ii[i]-1]=expos_time[jj[i]-100,ii[i]-1]+et[i]
    expos_part[jj[i]-100,ii[i]-1]=expos_part[jj[i]-100,ii[i]-1]+ep[i]

#-------------------------------------------------------------------------------
# Mask the age variables to the shelf
#-------------------------------------------------------------------------------
scope2=scope+hmask
scope=1-scope
hmask=1-hmask
cm.gnuplot2_r.set_bad('white',alpha=0)

#-------------------------------------------------------------------------------
# Open age file
#-------------------------------------------------------------------------------
day_count=[0,31,28,31,30,31,30,31,31,30,31,30,31]
fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor='w')
cm.gnuplot2_r.set_bad('w')

for i in range(0,1):

  #-----------------------------------------------------------------------------
  # Load and average data for each season
  #-----------------------------------------------------------------------------
  #print 'Plotting data for Month %2d' % (i+1)
  Ttime=expos_part*Gscope[99:260,0:200]
  Ttime=Ttime/(29230.0*6.0)

  #-----------------------------------------------------------------------------
  # Mask the age variables to the shelf
  #-----------------------------------------------------------------------------
  Ttime=np.ma.array(Ttime,mask=scope[99:260,0:200])
  Ttime=np.ma.masked_where(expos_part==0,Ttime) #<-- Mask where less than 0.5%
  
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
  cs = m.pcolormesh(xi,yi,Ttime,cmap=cm.gnuplot2_r,vmin=0, vmax = 0.15)
  xi, yi = m(lon,lat)
  cM = m.contour(   xi,yi,scope2,[1,1],colors='k',linewidth=2.0)
  #  xi, yi = m(lon[99:260,0:200],lat[99:260,0:200])
  #cn = m.contour(   xi,yi,Ttime,[182.625,365.25,730.5,913.125,1095.75],colors='0.5',linewidth=1.0)

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
  clabels=np.array([0,0.025,0.05,0.075,0.1,0.125,0.15])
  cbar = m.colorbar(cs, location='bottom', pad="4%",ticks=clabels,extend='max')
  cbar.ax.set_xticklabels([0,2.5,5,7.5,10,12.5,15])#,2.75,3])
  cbar.set_label('Probability of Particle Transit (%)')
  
  # Add Title
  #  plt.title('Season')

  # Save figure to file
  #plot_file=plot_dir+'month_'+str(i+1)+'.png'
  plot_file=plot_dir+'particle_PDF.png'
  print 'Saving figure to file %s ' % plot_file
  plt.savefig(plot_file,dpi=150,bbox_inches='tight')

