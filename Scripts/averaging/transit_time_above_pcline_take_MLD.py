#===============================================================================
#
# SCRIPT :  transit_time_above_pcline_take_MLD.py
#
# PURPOSE : Ingest input data from pycnocline-averaged US East water age model
#           and use the calcluated thermocline depth at each time to average the
#           residence time and velocities above the mixed layer depth.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 6 JULY 2015 for personal use
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
print " RUNNING: Script transit_time_above_pcline_take_MLD.py"
print "----------------------------------------------------------------------"
print ""

################################################################################
#                                                                              #
#                            User Defined Variables                            #
#                                                                              #
input_file  = '/Volumes/Backup2/Data/USEast-age/averages/avg_3hrly.nc'
MLD_file    = '/Volumes/Backup2/Data/USEast-age/averages/age_pycnocline_avg.nc'
Grtime_file = '/Volumes/Backup1/Data/USeast-rtime/rtime_GOM_shelf.nc'
Artime_file = '/Volumes/Backup1/Data/USeast-rtime/rtime_ATL.nc'
ATL_file    = '/Volumes/Backup2/Data/USEast/grid/grid_ATLscope.nc'
GOM_file    = '/Volumes/Backup2/Data/USEast/grid/grid_GOM_shelf_scope.nc'
grid_file   = '/Volumes/Backup2/Data/USEast/grid/USeast-grid.nc'
output_file = '/Volumes/Backup1/Data/USeast-rtime/rtime_pycnocline.nc'
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Open grid file and read in variables
#-------------------------------------------------------------------------------
print 'Ingesting data from file %s ...' % (grid_file)
grid_data=Dataset(grid_file,mode='r')
mask=grid_data.variables['mask_rho' ][:,:]
lon =grid_data.variables['lon_rho'  ][:,:]
lat =grid_data.variables['lat_rho'  ][:,:]
h   =grid_data.variables['h'        ][:,:]
grid_data.close()

print 'Ingesting data from file %s ...' % (ATL_file)
grid_data=Dataset(ATL_file,mode='r')
hmask  = grid_data.variables['mask_rho' ][:]
Ascope = grid_data.variables['scope_rho'][:]
grid_data.close()

# Get GOM Shelf gird
print 'Ingesting data from file %s ...' % (GOM_file)
grid_data=Dataset(GOM_file,mode='r')
Gscope = grid_data.variables['scope_rho'][:]
grid_data.close()

scope=Ascope+Gscope

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

hc=200.

#-------------------------------------------------------------------------------
# Open residence time files
#-------------------------------------------------------------------------------
Atime_data =Dataset(Artime_file,mode='r')
Gtime_data =Dataset(Grtime_file,mode='r')

#-------------------------------------------------------------------------------
# Open pycnocline averaged files
#-------------------------------------------------------------------------------
MLD_data =Dataset(MLD_file,mode='r')

#--------------------------------------------------------------------------
#  Create local storage variable
#--------------------------------------------------------------------------
rtime01 = np.zeros(shape=(len(input_data.variables['ocean_time'])-1,M,L))

#--------------------------------------------------------------------------
# Loop through the water age output to calculate thermocline depth
#--------------------------------------------------------------------------
print "-----------------------------"
print "Beginning Main Loop Iteration"
print "-----------------------------"
for tdim in range(0,len(input_data.variables['ocean_time'])-1):

  print 'Time : %04i' % tdim

  #------------------------------------------------------------------------
  # Load in the SSH, Temp, and Salt to calculate layer & thermocline depths
  #------------------------------------------------------------------------
  zeta=input_data.variables['zeta'       ][tdim,:,:  ]

  #------------------------------------------------------------------------
  # Calculate depths of each rho layer
  #------------------------------------------------------------------------
  z  = np.zeros(shape=(N  ,M,L))
  zw = np.zeros(shape=(N+1,M,L))
  z  = depths_rho(h,zeta,s_r,Cs_r,hc)
  zw = depths_w(  h,zeta,s_w,Cs_w,hc)
  Hz = abs(zw[1:N+1,:,:]-zw[0:N,:,:])

  #------------------------------------------------------------------------
  # Find layer with max dT/d(sigma)
  #------------------------------------------------------------------------
  for i in range(0,M):
    for j in range(0,L):
      
      # Cell is over the water
      if scope[i,j]>0:
        
        
        #-----------------------------------------------------------------
        # Load in residence times and MLD
        #----------------------------------------------------------------
        Grtime=Gtime_data.variables['age_02'][tdim,:,i,j]
        Artime=Atime_data.variables['age_02'][tdim,:,i,j]
        MLD   =  MLD_data.variables['MLD'   ][tdim,i,j]
        rtime =Grtime+Artime

        #-----------------------------------------------------------------
        # Read depth of the pycnocline
        #----------------------------------------------------------------
        drdZ_ind=np.argmin(abs(abs(z[:,i,j])-MLD))
        
        #-----------------------------------------------------------------
        # Average age above the pycnocline
        #----------------------------------------------------------------
        rtime01[tdim,i,j]=sum(np.squeeze(   Hz[drdZ_ind::,i,j])*\
                                         rtime[drdZ_ind::])/\
                              np.squeeze(abs(z[drdZ_ind,i,j]))
      # Cell is not within scope
      else:
        rtime01[tdim,i,j] = np.nan

#-------------------------------------------------------------------------------
# Open output file and create dimensions and variables
#-------------------------------------------------------------------------------
print " "
print "Creating output file: %s" % output_file
print " "

output_data=Dataset(output_file,mode='w',format='NETCDF4')
output_data.createDimension('ocean_time',len(input_data.dimensions['ocean_time'])-1)
output_data.createDimension('xi_rho'    ,len(input_data.dimensions['xi_rho'    ]))
output_data.createDimension('eta_rho'   ,len(input_data.dimensions['eta_rho'   ]))

time_out= output_data.createVariable('ocean_time','float32',('ocean_time'),\
                                     zlib=True,complevel=9,shuffle=True)
lon_rho = output_data.createVariable('lon_rho','float32',('eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
lat_rho = output_data.createVariable('lat_rho','float32',('eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
rtime_1 = output_data.createVariable('GOM_shelf_rtime','float32',\
                                     ('ocean_time','eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)

#--------------------------------------------------------------------------
#  Write output data to file
#--------------------------------------------------------------------------
lon_rho[ :] = lon[:]
lat_rho[ :] = lat[:]
time_out[:] = input_data.variables['ocean_time'][0:2922]
rtime_1[:]  = rtime01

# Close Files
print 'Closing files'
input_data.close()
Atime_data.close()
Gtime_data.close()
MLD_data.close()
output_data.close()




