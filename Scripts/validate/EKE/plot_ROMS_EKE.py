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

input_file = '/Volumes/Black_box/Data/USeast-age/output/clim/analysis/eke_sla.nc'
plot_dir   = '/Volumes/Black_box/Data/PLOTS/USeast-age/validate/'
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
hmask = grid_data.variables['mask_rho'][:]
pmask = grid_data.variables['mask_psi'][:]
lon   = grid_data.variables['lon_rho' ][:]
lat   = grid_data.variables['lat_rho' ][:]
lonp  = grid_data.variables['lon_psi' ][:]
latp  = grid_data.variables['lat_psi' ][:]
grid_data.close()

lon2=lonp-0.5*(lonp[1,2]-lonp[1,1])
lat2=latp-0.5*(latp[2,1]-latp[1,1])

#-------------------------------------------------------------------------------
# Load ROMS EKE data from file
#-------------------------------------------------------------------------------
print 'Ingesting data from file %s ...' % (input_file)
input_data=Dataset(input_file,mode='r')
eke  = input_data.variables['eke'      ][:]
ekep = input_data.variables['eke_prime'][:]
input_data.close()

#-------------------------------------------------------------------------------
# Mask the age variables to the shelf
#-------------------------------------------------------------------------------
hmask=1-hmask

#-------------------------------------------------------------------------------
# Plot EKE Prime
#-------------------------------------------------------------------------------
print 'Plotting EKEprime data.'

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
xi, yi = m(lon2,lat2)
cs = m.pcolormesh(xi,yi,np.log10(ekep),cmap=cm.gist_ncar,vmin=-2.5, vmax = 0)

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
plot_file=plot_dir+'ROMS_EKE.png'
print 'Saving figure to file %s ' % plot_file
plt.savefig(plot_file,dpi=150,bbox_inches='tight')


