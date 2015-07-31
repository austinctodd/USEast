#===============================================================================
#
# SCRIPT :  cross-shelf-fluxes.py
#
# PURPOSE : Ingest velocity data from US East model and calculate the cross-shelf
#           fluxes from different shelf domain regions (GOM,SAB,MAB,GOME). 
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 11 July 2015 for personal use.
#
#===============================================================================

ATL_file = '/Volumes/Black_box/Data/USeast/Data/grd/grid_ATLscope.nc'
GOM_file = '/Volumes/Black_box/Data/USeast/Data/grd/grid_GOM_shelf_scope.nc'
fwd_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
           'avg_3hrly.nc'
shelf_pts= 'shelf_points.txt'
output_file='/Volumes/Black_box/Data/USeast-age/analysis/cross_shelf_fluxes.nc'

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
import gsw

#-------------------------------------------------------------------------------
# Define all required sub-libraries, sub-routines, and sub-modules
#-------------------------------------------------------------------------------
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from matplotlib import cm

#-------------------------------------------------------------------------------
# Open grid file and read in variables
#-------------------------------------------------------------------------------
print 'Ingesting data from file %s ...' % (ATL_file)
grid_data=Dataset(ATL_file,mode='r')
hmask  = grid_data.variables['mask_rho' ][:,:]
Ascope = grid_data.variables['scope_rho'][:,:]
h      = grid_data.variables['h'        ][:,:]
lon    = grid_data.variables['lon_rho'  ][:,:]
lat    = grid_data.variables['lat_rho'  ][:,:]
grid_data.close()

print 'Ingesting data from file %s ...' % (GOM_file)
grid_data=Dataset(GOM_file,mode='r')
Gscope = grid_data.variables['scope_rho'][:,:]
grid_data.close()

# Add scopes together to obtain shelf scope mask
Sscope=Ascope+Gscope

# Shift lat/lons by 1/2 distance for plotting purposes
lon2=lon-0.5*(lon[1,2]-lon[1,1])
lat2=lat-0.5*(lat[2,1]-lat[1,1])

#-------------------------------------------------------------------------------
# Open forward file and read in static variables
#-------------------------------------------------------------------------------
fwd_data=Dataset(fwd_file,'r')

s_r =fwd_data.variables['s_rho'][:]
s_w =fwd_data.variables['s_w'  ][:]
Cs_r=fwd_data.variables['Cs_r' ][:]
Cs_w=fwd_data.variables['Cs_w' ][:]

L=len(fwd_data.dimensions['xi_rho' ])
M=len(fwd_data.dimensions['eta_rho'])
N=len(fwd_data.dimensions['s_rho'  ])

hc=200.

#-------------------------------------------------------------------------------
# Read in Shelf edge index points from file (skip first line)
#-------------------------------------------------------------------------------
shelf_i,shelf_j=np.loadtxt(shelf_pts,unpack=True,skiprows=1)
shelf_i=shelf_i.astype(np.int64)
shelf_j=shelf_j.astype(np.int64)

uflx=np.empty((2923,N,len(shelf_i),))
vflx=np.empty((2923,N,len(shelf_i),))
uflx[:]=np.nan
vflx[:]=np.nan

#-------------------------------------------------------------------------------
# Loop through each point calculate the angle of the shelf break at each point
#-------------------------------------------------------------------------------
for i in range (0,len(shelf_i)):
  
    print 'Fluxes for point %03i' % i
    
    #------------------------------------------------------------------------
    # Read in velocity and sea level data
    #------------------------------------------------------------------------
    gsw.distance(lon[shelf_j[i],shelf_i[i])
    
    #------------------------------------------------------------------------
    # Calculate depths of each rho layer
    #------------------------------------------------------------------------
    z =np.zeros(shape=(2923,N  ,))
    zw=np.zeros(shape=(2923,N+1,))
    for k in range(0,len(s_r)):
        z0  =(hc*s_r[k]+Cs_r[k]*h[shelf_j[i],shelf_i[i]])/(hc+h[shelf_j[i],shelf_i[i]]);
        z[:,k]=zeta+(zeta+h[shelf_j[i],shelf_i[i]])*z0;
    
        z0   =(hc*s_w[k]+Cs_w[k]*h[shelf_j[i],shelf_i[i]])/(hc+h[shelf_j[i],shelf_i[i]]);
        zw[:,k]=zeta+(zeta+h[shelf_j[i],shelf_i[i]])*z0;

    # Add last depth for zw
    z0   =(hc*s_w[N]+Cs_w[N]*h[shelf_j[i],shelf_i[i]])/(hc+h[shelf_j[i],shelf_i[i]]);
    zw[:,N]=zeta+(zeta+h[shelf_j[i],shelf_i[i]])*z0;
    Hz = abs(zw[:,1:N+1]-zw[:,0:N])
    
    #------------------------------------------------------------------------
    # Multiply velocities by layer depths
    #------------------------------------------------------------------------
    uflx[:,:,i]=u*Hz
    vflx[:,:,i]=v*Hz

#-------------------------------------------------------------------------------
# Write data to file
#-------------------------------------------------------------------------------
print " "
print "Creating output file: %s" % output_file
print " "

output_data=Dataset(output_file,mode='w',format='NETCDF4')
output_data.createDimension('ocean_time',uflx.shape[0])
output_data.createDimension('s_rho'   ,len(s_r))
output_data.createDimension('shelf_pts',len(shelf_i))

u_flx  = output_data.createVariable('U_flux','float32',('ocean_time','s_rho','shelf_pts'),\
                                     zlib=True,complevel=9,shuffle=True)
v_flx  = output_data.createVariable('V_flux','float32',('ocean_time','s_rho','shelf_pts'),\
                                    zlib=True,complevel=9,shuffle=True)

u_flx[:] = uflx
v_flx[:] = vflx

output_data.close()