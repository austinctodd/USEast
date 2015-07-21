#===============================================================================
#
# SCRIPT :  avg_thermocline_age.py
#
# PURPOSE : Ingest input data from US East water age model and calculate the
#           thermocline depth at each time. Then, average the water age above
#           and below the depth.  Or should I just output the depths/indexes?
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
print " RUNNING: Script transit_time_above_pcline.py"
print "----------------------------------------------------------------------"
print ""

################################################################################
#                                                                              #
#                            User Defined Variables                            #
#                                                                              #
input_file  = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
              'avg_3hrly.nc'
Grtime_file = '/Volumes/Black_box/Data/USeast-rtime/output/GOM_shelf/'+\
              'rtime_GOM_shelf.nc'
Artime_file = '/Volumes/Black_box/Data/USeast-rtime/output/ATL/'+\
              'rtime_ATL.nc'
ATL_file    = '/Volumes/Black_box/Data/USeast/Data/grd/grid_ATLscope.nc'
GOM_file    = '/Volumes/Black_box/Data/USeast/Data/grd/grid_GOM_shelf_scope.nc'
grid_file   = '/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc'
output_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
              'pcline_avgs.nc'
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

#--------------------------------------------------------------------------
#  Create local storage variable
#--------------------------------------------------------------------------
rtime02 = np.zeros(shape=(len(input_data.variables['ocean_time']),M,L))

#--------------------------------------------------------------------------
# Loop through the water age output to calculate thermocline depth
#--------------------------------------------------------------------------
print "-----------------------------"
print "Beginning Main Loop Iteration"
print "-----------------------------"
for tdim in range(0,len(input_data.variables['ocean_time'])):

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
  #z  = z[::-1,:,:]
  zw = depths_w(  h,zeta,s_w,Cs_w,hc)
  #z  = z[::-1,:,:]
  Hz = abs(zw[1:N+1,:,:]-zw[0:N,:,:])

  #------------------------------------------------------------------------
  # Calculate various seawater properties
  #------------------------------------------------------------------------
  pres=sw.pres(z*-1,lat)

  #------------------------------------------------------------------------
  # Find layer with max dT/d(sigma)
  #------------------------------------------------------------------------
  dT=np.zeros(shape=(N-1))
  for i in range(0,M):
    for j in range(0,L):
      
      # Cell is over the water
      if mask[i,j]>0:
        #-----------------------------------------------------------------
        # Load in temp,salt, age, residence time
        #----------------------------------------------------------------
        temp  =input_data.variables['temp'       ][tdim,:,i,j]
        salt  =input_data.variables['salt'       ][tdim,:,i,j]
        Grtime=Gtime_data.variables['age_02'     ][tdim,:,i,j]
        Artime=Atime_data.variables['age_02'     ][tdim,:,i,j]
        rtime =Grtime+Artime

        #-----------------------------------------------------------------
        # Calculate depth of the pycnocline
        #----------------------------------------------------------------
        dens=sw.dens(salt,temp,np.squeeze(pres[:,i,j]))-1000.
        drdZ = abs(dens[1:N]-dens[0:N-1])/abs(np.squeeze(z[1:N,i,j]-z[0:N-1,i,j]))
        drdZ_ind=np.argmax(drdZ)
        
        #-----------------------------------------------------------------
        # Average residence time above the pycnocline
        #----------------------------------------------------------------
        # Cell is within shelf scope
        if Gscope[i,j]>0:
          rtime02[tdim,i,j]=sum(np.squeeze(   Hz[drdZ_ind::,i,j])*\
                                           rtime[drdZ_ind::])/\
                                np.squeeze(abs(z[drdZ_ind,i,j]))
        else:
          rtime02[tdim,i,j]=np.nan
    
      # Cell is over land
      else:
        rtime02[tdim,i,j] = np.nan

#-------------------------------------------------------------------------------
# Open output file and create dimensions and variables
#-------------------------------------------------------------------------------
print " "
print "Creating output file: %s" % output_file
print " "

output_data=Dataset(output_file,mode='a')
rtime_2 = output_data.createVariable('full_rtime','float32',\
                                     ('ocean_time','eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)

#--------------------------------------------------------------------------
#  Write output data to file
#--------------------------------------------------------------------------
rtime_2[:]  = rtime02

# Close Files
print 'Closing files'
input_data.close()
Atime_data.close()
Gtime_data.close()
output_data.close()




