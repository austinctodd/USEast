#===============================================================================
#
# SCRIPT :  plot_EKE.py
#
# PURPOSE : Plot maps of EKE from US East water age model and AVISO SLA data.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 27 MARCH 2015.
#
#===============================================================================

input_dir = '/Users/actodd/Documents/ROMS-utils/USeast-age/validate/EKE/'
plot_dir  = '/Volumes/Black_box/Data/PLOTS/USeast-age/validate/'
GOM_file   = '/Volumes/Black_box/Data/USeast/Data/grd/grid_GOM_shelf_scope.nc'

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script plot_EKE.py"
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
lon   = grid_data.variables['lon_rho' ][:]
lat   = grid_data.variables['lat_rho' ][:]
grid_data.close()

#-------------------------------------------------------------------------------
# Load AVISO EKE data from file
#-------------------------------------------------------------------------------
file=input_dir+'AVISO_eke_sla.txt'
Alon,Alat,AEKE=np.loadtxt(file,unpack=True)

AVISO_EKE=np.zeros(shape=(161,165))
AVISO_lon=np.zeros(shape=(161,165))
AVISO_lat=np.zeros(shape=(161,165))

count=0
for i in range(0,161):
  for j in range(0,165):
    AVISO_EKE[i,j]=AEKE[count]
    AVISO_lon[i,j]=Alon[count]
    AVISO_lat[i,j]=Alat[count]
    count=count+1

#-------------------------------------------------------------------------------
# Mask the age variables to the shelf
#-------------------------------------------------------------------------------
cm.hot.set_bad('w')

#-------------------------------------------------------------------------------
# Plot AVISO PDF
#-------------------------------------------------------------------------------
print 'Plotting AVISO data.'

fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor='w')
plt.clf()
m=Basemap(projection='merc',llcrnrlat=lat[0,0],urcrnrlat=lat[481,401-25],\
          llcrnrlon=lon[0,0],urcrnrlon=lon[481,401-25],lat_ts=20,resolution='f')

# Draw a think border around the whole map
m.drawmapboundary(linewidth=3)

#-----------------------------------------------------------------------------
# Mask the age variables to the shelf
#-----------------------------------------------------------------------------
#AVISO_EKE=np.ma.array(AVISO_EKE)
#AVISO_EKE=np.ma.masked_where(AVISO_EKE==0,AVISO_EKE) #<-- Mask where less than 0.5%

# Plot Data
lon2=AVISO_lon-0.5*(AVISO_lon[1,2]-AVISO_lon[1,1])
lat2=AVISO_lat-0.5*(AVISO_lat[2,1]-AVISO_lat[1,1])
xi, yi = m(lon2,lat2)
cs = m.pcolormesh(xi,yi,np.log10(AVISO_EKE),cmap=cm.gist_ncar,vmin=-2.5, vmax = 0)

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
clabels=np.array([-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0])
cbar = m.colorbar(cs, location='bottom', pad="4%",ticks=clabels,extend='both')
#cbar.ax.set_xticklabels([2,4,6,8,10])#,2.75,3])
#cbar.ax.set_xticklabels([0,2.5,5,7.5,10,12.5,15])#,2.75,3])
cbar.set_label(r'$\log$'+'(EKE) (m'+r'$^2$'+'/s'+r'$^2$'+')')

# Save figure to file
plot_file=plot_dir+'AVISO_EKE_sla.png'
print 'Saving figure to file %s ' % plot_file
plt.savefig(plot_file,dpi=150,bbox_inches='tight')


