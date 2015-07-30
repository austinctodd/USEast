#===============================================================================
#
# SCRIPT :  ROMS_EKE_sla.py
#
# PURPOSE : Ingest zeta from US East water age forward model and average over
#           the six years of forward model state.  
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 2 April 2015 for personal use. Modified
#           23 July 2015.
#
#===============================================================================

#-------------------------------------------------------------------------------
# Define all required libraries, routines, and modules
#-------------------------------------------------------------------------------
import numpy as np
import scipy.signal as signal
import netCDF4
import sys
import os

#-------------------------------------------------------------------------------
# Define all required sub-libraries, sub-routines, and sub-modules
#-------------------------------------------------------------------------------
from netCDF4 import Dataset
from math import pi

#-------------------------------------------------------------------------------
# Print message to screen
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script ROMS_EKE_sla.py"
print "----------------------------------------------------------------------"
print ""

################################################################################
#                                                                              #
#                            User Defined Variables                            #
#                                                                              #
infile1 = '/Volumes/Backup2/Data/USeast-age/fwdrun/part1/useast_his.nc'
infile2 = '/Volumes/Backup2/Data/USeast-age/fwdrun/part2/useast_his.nc'
infile3 = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
          'mean_vels.nc'
outfile = '/Volumes/Black_box/Data/USeast-age/output/clim/analysis/'+\
          'eke_sla_filtered100hrs.nc'
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Open input and output datasets
#-------------------------------------------------------------------------------
indata1=Dataset(infile1,mode='r')
indata2=Dataset(infile2,mode='r')
indata3=Dataset(infile3,mode='r')
outdata=Dataset(outfile,mode='w')

#-------------------------------------------------------------------------------
# Read grid variables
#-------------------------------------------------------------------------------
print 'Ingesting grid data from file %s ...' % (infile1)
mask = indata1.variables['mask_rho'][:]
lon   = indata1.variables['lon_rho' ][:]
lat   = indata1.variables['lat_rho' ][:]

#-------------------------------------------------------------------------------
# Define output dimensions
#-------------------------------------------------------------------------------
print 'Defining variable dimensions ...'
dims=('xi_psi','eta_psi')
for dim in dims:
  outdata.createDimension(dim,len(indata1.dimensions[dim]))
  
#-------------------------------------------------------------------------------
# Create variables
#-------------------------------------------------------------------------------
print 'Creating eke variable ...'
ekeout = outdata.createVariable('eke',indata1.variables['u'].dtype,\
                             ('eta_psi','xi_psi'))
ekeout.long_name = 'Mean EKE calculated from deviations from six year mean SSH'+\
                   '- filtered with 100-hr lowpass filter'
ekeout.units = 'm^2/s^2'
ekeout.field = 'mean eddy kinetic energy, scalar, series'

################################################################################
##                                                                            ##
##                         Average through velocities                         ##
##                                                                            ##
################################################################################

#--------------------------------------------------------------------------
# Declare constants used for calculations
#--------------------------------------------------------------------------
DEG2RAD = (2*pi/360.0);
DEG2NM  = 60.0;
NM2KM   = 1.8520;    # Defined in Pond & Pickard p303.
g       = -9.81;

#--------------------------------------------------------------------------
# Calculate meridional distances (dy)
#--------------------------------------------------------------------------
latrad = (indata1.variables['lat_rho'][:,:])*DEG2RAD
dlon   =  indata1.variables['lon_rho'][1:482,:]-\
          indata1.variables['lon_rho'][0:481,:]
dlat   =  indata1.variables['lat_rho'][1:482,:]-\
          indata1.variables['lat_rho'][0:481,:]
dep    = np.cos((latrad[1:482,:]+latrad[0:481,:])/2.0 )*dlon;
dy     = np.sqrt(dep**2 + dlat**2)*DEG2NM*NM2KM*1000.0
fy     = (indata1.variables['f'][1:482,:]+\
          indata1.variables['f'][0:481,:])/2.0

#--------------------------------------------------------------------------
# Calculate zonal distances (dx)
#--------------------------------------------------------------------------
latrad = (indata1.variables['lat_rho'][:,:])*DEG2RAD
dlon   =  indata1.variables['lon_rho'][:,1:402]-\
          indata1.variables['lon_rho'][:,0:401]
dlat   =  indata1.variables['lat_rho'][:,1:402]-\
          indata1.variables['lat_rho'][:,0:401]
dep    = np.cos((latrad[:,1:402]+latrad[:,0:401])/2.0 )*dlon;
dx     = np.sqrt(dep**2 + dlat**2)*DEG2NM*NM2KM*1000.0
fx     = (indata1.variables['f'][:,1:402]+\
          indata1.variables['f'][:,0:401])/2.0

#--------------------------------------------------------------------------
# Declare 2d variables to be averaged
#--------------------------------------------------------------------------
M=len(indata1.dimensions[indata1.variables['zeta'].dimensions[1]])
N=len(indata1.dimensions[indata1.variables['zeta'].dimensions[2]])
eke       = np.zeros(shape=(M-1,N-1))

#--------------------------------------------------------------------------
# 6th Order Butterworth filter w cutoff frequency of 100 hrs @ 3 hrly data
#--------------------------------------------------------------------------
B,A=signal.butter(6,3.0/50.0,output='ba');

#------------------------------------------------------------------------
# Initialize some arrays
#------------------------------------------------------------------------
count=0
zeta1=np.zeros(shape=(2922*6,)) # Lower left  corner [i  ,j  ]
zeta2=np.zeros(shape=(2922*6,)) # Lower right corner [i  ,j+1]
zeta3=np.zeros(shape=(2922*6,)) # Upper left  corner [i+1,j  ]
zeta4=np.zeros(shape=(2922*6,)) # Upper right corner [i+1,j+1]

for i in range(0,lon.shape[0]-1):
  for j in range(0,lon.shape[1]-1):

    print 'Cell %i / %i' % (count,401*481)
    count+=1

    if (mask[i,j] > 0 and mask[i+1,j]>0 and mask[i,j+1]>0 and mask[i+1,j+1]>0):
      
      #----------------------------------------------------------------------
      # Calculate geostrophic velocities based on deviations from ROMS zeta:
      #                  z  = zbar - z (zbar is long-term zeta mean)
      #                  ug = -(g/f)*(dz/dy)
      #                  vg =  (g/f)*(dz/dx)
      #----------------------------------------------------------------------

      #----------------------------------------------------------------------
      # Load zeta at i,j
      #----------------------------------------------------------------------
      zeta1[0:2922*2]       =indata1.variables['zeta'][0:2922*2  ,i,j]
      zeta1[2922*2-1:2922*6]=indata2.variables['zeta'][0:2922*4+1,i,j]

      # Filter the SSH over the entire record
      zeta1f=signal.filtfilt(B,A,zeta1-indata3.variables['mean_zeta'][i,j])
      cff1=zeta1f[4:2922*6-4]

      #----------------------------------------------------------------------
      # Load zeta at i,j+1
      #----------------------------------------------------------------------
      zeta2[0:2922*2]       =indata1.variables['zeta'][0:2922*2  ,i,j+1]
      zeta2[2922*2-1:2922*6]=indata2.variables['zeta'][0:2922*4+1,i,j+1]

      # Filter the SSH over the entire record
      zeta2f=signal.filtfilt(B,A,zeta2-indata3.variables['mean_zeta'][i,j+1])
      cff2=zeta2f[4:2922*6-4]

      #----------------------------------------------------------------------
      # Load the variables at i+1,j
      #----------------------------------------------------------------------
      zeta3[0:2922*2]       =indata1.variables['zeta'][0:2922*2  ,i+1,j]
      zeta3[2922*2-1:2922*6]=indata2.variables['zeta'][0:2922*4+1,i+1,j]

      # Filter the SSH over the entire record
      zeta3f=signal.filtfilt(B,A,zeta3-indata3.variables['mean_zeta'][i+1,j])
      cff3=zeta3f[4:2922*6-4]

      #----------------------------------------------------------------------
      # Load the variables at i+1,j+1
      #----------------------------------------------------------------------
      zeta4[0:2922*2]       =indata1.variables['zeta'][0:2922*2  ,i+1,j+1]
      zeta4[2922*2-1:2922*6]=indata2.variables['zeta'][0:2922*4+1,i+1,j+1]
    
      # Filter the SSH over the entire record
      zeta4f=signal.filtfilt(B,A,zeta4-indata3.variables['mean_zeta'][i+1,j+1])
      cff4=zeta4f[4:2922*6-4]

      #----------------------------------------------------------------------
      # Calculate ug and vg at psi-points
      #----------------------------------------------------------------------
      ug = -(g/2.0)*((cff3-cff1)/(fy[i,j]*dy[i,j])+(cff4-cff2)/(fy[i,j+1]*dy[i,j+1]))
      vg=   (g/2.0)*((cff2-cff1)/(fx[i,j]*dx[i,j])+(cff4-cff3)/(fx[i+1,j]*dx[i+1,j]))

      # Now calculate the mean over each day
      eke[i,j]=np.nanmean(0.5*((ug**2) + (vg**2)))
    else:
      eke[i,j]=np.nan

# Output data
ekeout[:] = eke

indata1.close()
indata2.close()
indata3.close()
outdata.close()
