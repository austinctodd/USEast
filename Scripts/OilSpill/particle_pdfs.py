#===============================================================================
#
# SCRIPT :  particle_pdfs.py
#
# PURPOSE : Ingest data from LTRANS run of simulated Oil Spill release locations
#           and calculcate the probability density function that particles will
#           enter each grid cell of model domain.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 30 July 2015 for personal use.
#
#===============================================================================

grid_file  = '/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc'
LTRANS_dir = '/Volumes/Black_box/Data/LTRANS/output/alvinocaris/'
cmap_file  = '/Users/actodd/Documents/ROMS-utils/USeast/'+\
             'ElevationBathymetryColormap/etopo1.clr'
#plot_dir   = '/Users/todd/Documents/Work/Projects/Atlantis/'
output_file= '/Volumes/Black_box/Data/LTRANS/output/OilSpill/pcounts_14.nc'

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script particle_pdfs.py"
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
from matplotlib.patches import Polygon
#from numpy import append, arange, dtype, linspace, meshgrid, ndarray

#-------------------------------------------------------------------------------
# Open grid file and read in variables
#-------------------------------------------------------------------------------
print 'Ingesting data from file %s ...' % (grid_file)
grid_data=Dataset(grid_file,mode='r')
hmask = grid_data.variables['mask_rho'][:]
hraw = grid_data.variables['hraw'     ][0,:]
h     = grid_data.variables['h'       ][:]
lon   = grid_data.variables['lon_rho' ][:,:]
lat   = grid_data.variables['lat_rho' ][:,:]
grid_data.close()

# Make shifted lat & lon grids for plotting purposes
lon2=lon-0.5*(lon[1,2]-lon[1,1])
lat2=lat-0.5*(lat[2,1]-lat[1,1])

#-------------------------------------------------------------------------------
# Mask the variables over land
#-------------------------------------------------------------------------------
hraw=np.where(hraw>0.0,0.0,1.0)
for i in range(0,482):
  for j in range(0,401-25):
    if (hraw[i,j]==1.0 and hmask[i,j]==1.0):
        hraw[i,j]=0.0

h=np.ma.array(h,mask=hraw)

#-------------------------------------------------------------------------------
# Open each LTRANS file
#-------------------------------------------------------------------------------
particle_data=Dataset(LTRANS_dir+'site14.nc',mode='r')

#-----------------------------------------------------------------------------
# Loop particle data
#-----------------------------------------------------------------------------
pbox=np.zeros(shape=lon.shape)
PLD=365.0
for t in range(0,2925,4):
  stind=np.max([(t/4.0-PLD)*286.0,     0])
  enind=np.min([ t/4.0     *286.0,104676])
  
  if (stind >= 29280):
    break;
  else:
    plon=particle_data.variables['lon'  ][t,stind:enind]
    plat=particle_data.variables['lat'  ][t,stind:enind]
    pcol=particle_data.variables['color'][t,stind:enind]
    
    bb=np.where(np.array(pcol)<0)[0]
    plon[bb]=np.nan
    plat[bb]=np.nan
    
    for i in range(1,lon.shape[0]):
      for j in range(1,lon.shape[1]):
        if (hmask[i,j]>0 and hmask[i-1,j-1]>0 and hmask[i-1,j]>0 and hmask[i,j-1]>0):
          aa=plat[np.logical_and(np.array(plon)>=lon[i-1,j-1],np.array(plon)<lon[i,j])]
          bb=plat[np.logical_and(np.array(aa  )>=lat[i-1,j-1],np.array(aa  )<lat[i,j])]
          pbox[i,j]=pbox[i,j]+len(bb)

particle_data.close()

#-----------------------------------------------------------------------------
# Output data to file
#-----------------------------------------------------------------------------
print " "
print "Creating output file: %s" % output_file
print " "

output_data=Dataset(output_file,mode='w',format='NETCDF4')
output_data.createDimension('xi_rho'    ,lon.shape[1])
output_data.createDimension('eta_rho'   ,lon.shape[0])

lon_rho = output_data.createVariable('lon_rho','float32',('eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
lat_rho = output_data.createVariable('lat_rho','float32',('eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
p_box   = output_data.createVariable('counts','float32',('eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)

p_box[:] = pbox

output_data.close()

site14.close()
#site15.close()
#site16.close()
