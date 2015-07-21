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
#--------------------------------------------------------------------------
def depths_rho(h,zeta,s_r,Cs_r,hc):
  for k in range(0,len(s_r)):
    z0=(hc*s_r[k]+Cs_r[k]*h)/(hc+h);
    z[k,:,:]=zeta+(zeta+h)*z0;
  return z

#--------------------------------------------------------------------------
# Routine for calculating depths at u points
#--------------------------------------------------------------------------
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
# Routine for calculating depths at v points
#--------------------------------------------------------------------------
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
# Routine for calculating depths at w points
#--------------------------------------------------------------------------
def depths_w(h,zeta,s_w,Cs_w,hc):
  for k in range(0,len(s_w)):
    z0=(hc*s_w[k]+Cs_w[k]*h)/(hc+h);
    zw[k,:,:]=zeta+(zeta+h)*z0;
  return zw

#-------------------------------------------------------------------------------
# Print message to screen
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script age_above_pcline.py"
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

#-------------------------------------------------------------------------------
# Open output file and create dimensions and variables
#-------------------------------------------------------------------------------
print " "
print "Creating output file: %s" % output_file
print " "

output_data=Dataset(output_file,mode='w',format='NETCDF4')
output_data.createDimension('ocean_time',len(input_data.dimensions['ocean_time']))
output_data.createDimension('xi_rho'    ,len(input_data.dimensions['xi_rho'    ]))
output_data.createDimension('xi_u'      ,len(input_data.dimensions['xi_u'      ]))
output_data.createDimension('xi_v'      ,len(input_data.dimensions['xi_v'      ]))
output_data.createDimension('eta_rho'   ,len(input_data.dimensions['eta_rho'   ]))
output_data.createDimension('eta_u'     ,len(input_data.dimensions['eta_u'     ]))
output_data.createDimension('eta_v'     ,len(input_data.dimensions['eta_v'     ]))

time_out= output_data.createVariable('ocean_time','float32',('ocean_time'),\
                                     zlib=True,complevel=9,shuffle=True)
lon_rho = output_data.createVariable('lon_rho','float32',('eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
lat_rho = output_data.createVariable('lat_rho','float32',('eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
age_01  = output_data.createVariable('mean_age_01','float32',\
                                     ('ocean_time','eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
age_02  = output_data.createVariable('mean_age_02','float32',\
                                     ('ocean_time','eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
age_03  = output_data.createVariable('mean_age_03','float32',\
                                     ('ocean_time','eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
mld     = output_data.createVariable('MLD','float32',\
                                     ('ocean_time','eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
#u_out   = output_data.createVariable('u','float32',\
#                                     ('ocean_time','s_rho','eta_u','xi_u'),\
#                                     zlib=True,complevel=9,shuffle=True)
#v_out   = output_data.createVariable('v','float32',\
#                                     ('ocean_time','s_rho','eta_v','xi_v'),\
#                                     zlib=True,complevel=9,shuffle=True)
mixld = np.zeros(shape=(len(input_data.variables['ocean_time']),M,L))
age01 = np.zeros(shape=(len(input_data.variables['ocean_time']),M,L))
age02 = np.zeros(shape=(len(input_data.variables['ocean_time']),M,L))
age03 = np.zeros(shape=(len(input_data.variables['ocean_time']),M,L))
                 
#--------------------------------------------------------------------------
#  Create standard variables (lon,lat,etc.)
#--------------------------------------------------------------------------
lon_rho[ :] = input_data.variables['lon_rho'   ][:]
lat_rho[ :] = input_data.variables['lat_rho'   ][:]
time_out[:] = input_data.variables['ocean_time'][:]

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
  temp=input_data.variables['temp'       ][tdim,:,:,:]
  salt=input_data.variables['salt'       ][tdim,:,:,:]
  age1=input_data.variables['mean_age_01'][tdim,:,:,:]
  age2=input_data.variables['mean_age_02'][tdim,:,:,:]
  age3=input_data.variables['mean_age_03'][tdim,:,:,:]
  u   =input_data.variables['u'          ][tdim,:,:,:]
  v   =input_data.variables['v'          ][tdim,:,:,:]

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
  dens=sw.dens(salt,temp,pres)-1000.

  #------------------------------------------------------------------------
  # Find layer with max dT/d(sigma)
  #------------------------------------------------------------------------
  dT=np.zeros(shape=(N-1))
  for i in range(0,M):
    for j in range(0,L):
      
      # Cell is over the water
      if mask[i,j]>0:
        drdZ = abs(dens[1:N-1,i,j]-dens[0:N-2,i,j])/abs(z[1:N-1,i,j]-z[0:N-2,i,j])
        drdZ_ind=np.argmax(drdZ)
        
        #-----------------------------------------------------------------
        # Average above and below the pycnocline
        #----------------------------------------------------------------
        mixld[tdim,i,j]=abs(   z[drdZ_ind,   i,j])
        age01[tdim,i,j]=sum(  Hz[drdZ_ind:35,i,j]*\
                            age1[drdZ_ind:35,i,j])/abs(z[drdZ_ind  ,i,j])
        age02[tdim,i,j]=sum(  Hz[drdZ_ind:35,i,j]*\
                            age2[drdZ_ind:35,i,j])/abs(z[drdZ_ind  ,i,j])
        age03[tdim,i,j]=sum(  Hz[drdZ_ind:35,i,j]*\
                            age3[drdZ_ind:35,i,j])/abs(z[drdZ_ind  ,i,j])

      # Cell is over land
      else:
        mixld[tdim,i,j] = np.nan
        age01[tdim,i,j] = np.nan
        age02[tdim,i,j] = np.nan
        age03[tdim,i,j] = np.nan

print "Writing to file"
mld[   :]=mixld
age_01[:]=age01
age_02[:]=age02
age_03[:]=age03

input_data.close()
output_data.close()




