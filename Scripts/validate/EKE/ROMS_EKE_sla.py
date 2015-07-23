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
          'eke_sla_filtered.nc'
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
hmask = indata1.variables['mask_rho'][:]
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
eke = outdata.createVariable('eke',indata1.variables['u'].dtype,\
                             ('eta_psi','xi_psi'))
eke.long_name = 'Mean EKE calculated from deviations from six year mean SSH'+\
                '- filtered with 100-hr lowpass filter'
eke.units = 'm^2/s^2'
eke.field = 'mean eddy kinetic energy, scalar, series'

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
masky  =  indata1.variables['mask_rho'][1:482,:]*\
          indata1.variables['mask_rho'][0:481,:]

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
maskx  =  indata1.variables['mask_rho'][:,1:402]*\
          indata1.variables['mask_rho'][:,0:401]

#--------------------------------------------------------------------------
# Declare 2d variables to be averaged
#--------------------------------------------------------------------------
M=len(indata1.dimensions[indata1.variables['zeta'].dimensions[1]])
N=len(indata1.dimensions[indata1.variables['zeta'].dimensions[2]])
eke       = np.zeros(shape=(M-1,N-1))

#--------------------------------------------------------------------------
# 6th Order Butterworth filter w cutoff frequency of 100 hrs @ 3 hrly data
#--------------------------------------------------------------------------
B,A=scipy.signal.butter(6,3/50,output='ba');

count=0
#------------------------------------------------------------------------
# Determine indexes corresponding to date from each file
#------------------------------------------------------------------------
inds=[0+tdim,2922+tdim,0+tdim,2922+tdim,5844+tdim,8766+tdim]

zeta=np.zeros(shape=(2922*6,))

for i in range(0,lon.shape[0]-2):
  for j in range(0,lon.shape[1]-2):

    if mask[i,j] > 0:
      
      #----------------------------------------------------------------------
      # Calculate geostrophic velocities based on deviations from ROMS zeta:
      #                  z  = zbar - z (zbar is long-term zeta mean)
      #                  ug = -(g/f)*(dz/dy)
      #                  vg =  (g/f)*(dz/dx)
      #----------------------------------------------------------------------

      #----------------------------------------------------------------------
      # Load the variables for i+1
      #----------------------------------------------------------------------
      zeta[0:2922*2]       =indata1.variables['zeta'][0:2922*2,i+1,j]
      zeta[2922*2-1:2922*6]=indata2.variables['zeta'][0:2922*4,i+1,j]

      # Filter the SSH over the entire record
      zetaf=signal.filtfilt(B,A,zeta)

      # Filter the SSH over the entire record
      cff1 = zeta-indata3.variables['mean_zeta'][i+1,j]

      #----------------------------------------------------------------------
      # Load the variables at j+1
      #----------------------------------------------------------------------
      zeta[0:2922*2]       =indata1.variables['zeta'][0:2922*2,i,j+1]
      zeta[2922*2-1:2922*6]=indata2.variables['zeta'][0:2922*4,i,j+1]

      # Filter the SSH over the entire record
      zetaf=signal.filtfilt(B,A,zeta)

      # Filter the SSH over the entire record
      cff2 = zeta-indata3.variables['mean_zeta'][i,j+1]
 
      #----------------------------------------------------------------------
      # Load the variables at i,j
      #----------------------------------------------------------------------
      zeta[0:2922*2]       =indata1.variables['zeta'][0:2922*2,i,j]
      zeta[2922*2-1:2922*6]=indata2.variables['zeta'][0:2922*4,i,j]

      # Filter the SSH over the entire record
      zetaf=signal.filtfilt(B,A,zeta)

      #----------------------------------------------------------------------
      # Finish ug and vg calculation
      #----------------------------------------------------------------------
      ug[i,j] = -(masky*g/fy/dy)*(cff1-(zetaf-indata3.variables['mean_zeta'][i,j]))
      vg[i,j]=   (maskx*g/fx/dx)*(cff2-(zetaf-indata3.variables['mean_zeta'][i,j]))
        
      eke[i,j]=eke_prime+0.5*(((ug[:,1:402]+ug[:,0:401])/2.0)**2 +\
                              ((vg[1:482,:]+vg[0:481,:])/2.0)**2)

  count=count+6
    
eke=eke/count
eke_prime=eke_prime/count

eke1[:]=eke
eke2[:]=eke_prime

indata1.close()
indata2.close()
indata3.close()
outdata.close()
