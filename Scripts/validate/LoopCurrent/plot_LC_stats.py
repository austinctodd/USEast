#===============================================================================
#
# SCRIPT :  plot_LC_stats.py
#
# PURPOSE : Ingest positions of the Loop Current from from US East water age
#           model and from AVISO SLA data. Plot PDF maps of LC position.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 27 MARCH 2015.
#
#===============================================================================

input_dir = '/Users/actodd/Documents/ROMS-utils/USeast-age/validate/LoopCurrent/'
plot_dir  = '/Volumes/Black_box/Data/PLOTS/USeast-age/validate/'
GOM_file   = '/Volumes/Black_box/Data/USeast/Data/grd/grid_GOM_shelf_scope.nc'

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script plot_LC_stats.py"
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
hmask = grid_data.variables['mask_rho' ][104:267,0:198]
scope = grid_data.variables['scope_rho'][104:267,0:198]
lon   = grid_data.variables['lon_rho'  ][104:267,0:198]
lat   = grid_data.variables['lat_rho'  ][104:267,0:198]
grid_data.close()

lon2=lon-0.5*(lon[1,2]-lon[1,1])
lat2=lat-0.5*(lat[2,1]-lat[1,1])

#-------------------------------------------------------------------------------
# Load Loop Current position data from file
#-------------------------------------------------------------------------------
file=input_dir+'AVISO_LC_positions.txt'
Aii,Ajj,ALC,Acnt=np.loadtxt(file,unpack=True)

file=input_dir+'ROMS_LC_positions.txt'
Rii,Rjj,RLC,Rcnt=np.loadtxt(file,unpack=True)

LC_AVISO=np.zeros(shape=hmask.shape)
LC_ROMS =np.zeros(shape=hmask.shape)

for i in range(0,Aii.size):
  LC_AVISO[Ajj[i]-1,Aii[i]-1]=LC_AVISO[Ajj[i]-1,Aii[i]-1]+ALC[i]
  LC_ROMS[ Rjj[i]-1,Rii[i]-1]=LC_ROMS[ Rjj[i]-1,Rii[i]-1]+RLC[i]

LC_AVISO=LC_AVISO/Acnt[0]
LC_ROMS =LC_ROMS /Rcnt[0]

#-------------------------------------------------------------------------------
# Mask the age variables to the shelf
#-------------------------------------------------------------------------------
scope2=scope+hmask
scope=1-scope
hmask=1-hmask
cm.gnuplot2_r.set_bad('w')

#-------------------------------------------------------------------------------
# Plot ROMS PDF
#-------------------------------------------------------------------------------
print 'Plotting ROMS data.'

fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor='w')
plt.clf()
m=Basemap(projection='merc',llcrnrlat=18,urcrnrlat=31,\
          llcrnrlon=-98.5,urcrnrlon=-80,lat_ts=20,resolution='f')

# Draw a think border around the whole map
m.drawmapboundary(linewidth=3)

# Plot Data
xi, yi = m(lon2,lat2)
cs = m.pcolormesh(xi,yi,LC_ROMS,cmap=cm.gnuplot2_r,vmin=0.01, vmax = 0.1)
xi, yi = m(lon,lat)
cM = m.contour(   xi,yi,scope2,[1,1],colors='k',linewidth=2.0)
  
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
#clabels=np.array([0,0.025,0.05,0.075,0.1,0.125,0.15])
clabels=np.array([0.02,0.04,0.06,0.08,0.1])
cbar = m.colorbar(cs, location='bottom', pad="4%",ticks=clabels,extend='max')
cbar.ax.set_xticklabels([2,4,6,8,10])#,2.75,3])
cbar.set_label('Probability of Loop Current Position (%)')
  
# Save figure to file
plot_file=plot_dir+'LC_ROMS.png'
print 'Saving figure to file %s ' % plot_file
plt.savefig(plot_file,dpi=150,bbox_inches='tight')


#-------------------------------------------------------------------------------
# Plot AVISO PDF
#-------------------------------------------------------------------------------
print 'Plotting AVISO data.'

fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor='w')
plt.clf()
m=Basemap(projection='merc',llcrnrlat=18,urcrnrlat=31,\
          llcrnrlon=-98.5,urcrnrlon=-80,lat_ts=20,resolution='f')

# Draw a think border around the whole map
m.drawmapboundary(linewidth=3)

# Plot Data
xi, yi = m(lon2,lat2)
cs = m.pcolormesh(xi,yi,LC_AVISO,cmap=cm.gnuplot2_r,vmin=0.01, vmax = 0.1)
xi, yi = m(lon,lat)
cM = m.contour(   xi,yi,scope2,[1,1],colors='k',linewidth=2.0)

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
#clabels=np.array([0,0.025,0.05,0.075,0.1,0.125,0.15])
clabels=np.array([0.02,0.04,0.06,0.08,0.1])
cbar = m.colorbar(cs, location='bottom', pad="4%",ticks=clabels,extend='max')
cbar.ax.set_xticklabels([2,4,6,8,10])#,2.75,3])
#cbar.ax.set_xticklabels([0,2.5,5,7.5,10,12.5,15])#,2.75,3])
cbar.set_label('Probability of Loop Current Position (%)')

# Save figure to file
plot_file=plot_dir+'LC_AVISO.png'
print 'Saving figure to file %s ' % plot_file
plt.savefig(plot_file,dpi=150,bbox_inches='tight')


