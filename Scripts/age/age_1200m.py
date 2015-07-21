#===============================================================================
#
# SCRIPT :  age_1200m.py
#
# PURPOSE : Ingest input data from US East water age model and calculate the
#           thermocline depth at each time. Then, average the water age above
#           and below the depth.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 26 JANUARY 2015.
#
#===============================================================================

#-------------------------------------------------------------------------------
# Define all required libraries, routines, and modules
#-------------------------------------------------------------------------------
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import sys
import mpl_toolkits.basemap
import os
#import pyroms
sys.path.append('/Users/actodd/MYPYTHON/seawater/')
import seawater as sw

#-------------------------------------------------------------------------------
# Define all required sub-libraries, sub-routines, and sub-modules
#-------------------------------------------------------------------------------
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
from scipy.interpolate import griddata
#from pyroms.depths import nc_depths

#--------------------------------------------------------------------------
# Routine for calculating depths at rho points
def depths_rho(h,zeta,s_r,Cs_r,hc):
  for k in range(0,len(s_r)):
    z0=(hc*s_r[k]+Cs_r[k]*h)/(hc+h);
    z[k,:,:]=zeta+(zeta+h)*z0;
  return z
#  endfor
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Routine for calculating depths at u points
def depths_u(h,zeta,s_r,Cs_r,hc):
  M=len(h)
  L=len(h[0])
  hu=   0.5*(   h[0:M-1,1:L-1]+   h[0:M-1,0:L-2]);
  zetau=0.5*(zeta[0:M-1,1:L-1]+zeta[0:M-1,0:L-2])
  
  # Loop through each depth layer
  for k in range(0,len(s_r)):
    z0=(hc*s_r[k]+Cs_r[k]*hu)/(hc+hu);
    z[k,:,:]=zetau+(zetau+hu)*z0;
  #  endfor
  
  del hu
  del zetau
  return z
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Routine for calculating depths at v points
def depths_v(h,zeta,s_r,Cs_r,hc):
  M=len(h)
  L=len(h[0])
  hv=   0.5*(   h[1:M-1,0:L-1]+   h[0:M-2,0:L-1]);
  zetav=0.5*(zeta[1:M-1,0:L-1]+zeta[0:M-2,0:L-1])
  
  # Loop through all layers
  for k in range(0,len(s_r)):
    z0=(hc*s_r[k]+Cs_r[k]*hv)/(hc+hv);
    z[k,:,:]=zetav+(zetav+hv)*z0;
  #  endfor
  del hv
  del zetav
  return z
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Routine for calculating depths at w points
def depths_w(h,zeta,s_w,Cs_w,hc):
  for k in range(0,len(s_w)):
    z0=(hc*s_w[k]+Cs_w[k]*h)/(hc+h);
    zw[k,:,:]=zeta+(zeta+h)*z0;
  return zw
#  endfor
#--------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Print message to screen
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script plot_age_panels.py"
print "----------------------------------------------------------------------"
print ""

################################################################################
#                                                                              #
#                            User Defined Variables                            #
#                                                                              #
input_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
             'avg_3hrly.nc'
grid_file  = '/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc'
output_file= '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
             'age_pycnocline_avg.nc'
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Open grid file and read in variables
#-------------------------------------------------------------------------------
grid_data=Dataset(grid_file,mode='r')
mask=grid_data.variables['mask_rho' ][:,:]
lon =grid_data.variables['lon_rho'  ][:,:]
lat =grid_data.variables['lat_rho'  ][:,:]
h   =grid_data.variables['h'        ][:,:]
grid_data.close()

hc=200.

#-------------------------------------------------------------------------------
# Read in additional grid variables from forward file
#-------------------------------------------------------------------------------
input_data=Dataset(input_file,mode='r')
s_r =input_data.variables['s_rho'][:]
s_w =input_data.variables['s_w'  ][:]
Cs_r=input_data.variables['Cs_r' ][:]
Cs_w=input_data.variables['Cs_w' ][:]

L=len(input_data.dimensions['xi_rho' ])
M=len(input_data.dimensions['eta_rho'])
N=len(input_data.dimensions['s_rho'  ])

#--------------------------------------------------------------------------
# Loop through the water age output to calculate thermocline depth
#--------------------------------------------------------------------------
print "-----------------------------"
print "Beginning Main Loop Iteration"
print "-----------------------------"
day_count=[0,31,28,31,30,31,30,31,31,30,31,30,31]

for i in range(0,1):

  print 'Season : %04i' % (i+1)

  #------------------------------------------------------------------------
  # Load in the SSH, Temp, and Salt to calculate layer & thermocline depths
  #------------------------------------------------------------------------
  stind = (sum(day_count[0:(i  )*3+1]))*8
  enind = (sum(day_count[0:(i+1)*3+1]))*8-1

  zeta=np.nanmean(input_data.variables['zeta'       ][stind:enind,:,:  ],axis=0)
  age2=np.nanmean(input_data.variables['mean_age_02'][stind:enind,:,:,:],axis=0)
#  age3=np.nanmean(input_data.variables['mean_age_03'][stind:enind,:,:,:],axis=0)

  #------------------------------------------------------------------------
  # Calculate depths of each rho layer
  #------------------------------------------------------------------------
  z  = np.zeros(shape=(N  ,M,L))
  z  = depths_rho(h,zeta,s_r,Cs_r,hc)
  zi = np.zeros(shape=(M,L))-1200.0

  #------------------------------------------------------------------------
  # Calculate various seawater properties
  #------------------------------------------------------------------------
  print 'interpolating'
  age_02=griddata(age2,z,zi,method='linear',fill_value=np.nan)

input_data.close()

print 'done'



