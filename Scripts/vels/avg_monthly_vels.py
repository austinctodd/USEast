#===============================================================================
#
# SCRIPT :  plot_mean_vels.py
#
# PURPOSE : Ingest input data from US East water age model and plot the monthly
#           mean velocities on a map.  Try masking out non-shelf values
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 4 MARCH 2015.
#
#===============================================================================

input_file  = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
              'avg_3hrly.nc'
output_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
              'monthly_avg_vels.nc'

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script avg_monthly_vels.py"
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
# Open input file
#-------------------------------------------------------------------------------
input_data=Dataset(input_file,mode='r')

#-------------------------------------------------------------------------------
# Open output file and create dimensions and variables
#-------------------------------------------------------------------------------
print " "
print "Creating output file: %s" % output_file
print " "

output_data=Dataset(output_file,mode='w',format='NETCDF4')
output_data.createDimension('ocean_time',12)
output_data.createDimension('s_rho'  ,36)
output_data.createDimension('xi_rho' ,400)
output_data.createDimension('eta_rho',480)

lon_rho = output_data.createVariable('lon_rho','float32',('eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
lat_rho = output_data.createVariable('lat_rho','float32',('eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)

#--------------------------------------------------------------------------
#  Create standard variables (lon,lat,etc.)
#--------------------------------------------------------------------------
lon_rho[:] = input_data.variables['lon_rho'][1:481,1:401]
lat_rho[:] = input_data.variables['lat_rho'][1:481,1:401]

#-------------------------------------------------------------------------------
# Define output variables
#-------------------------------------------------------------------------------
day_count=[0,31,28,31,30,31,30,31,31,30,31,30,31]

uu=np.zeros(shape=(12,36,480,400))
vv=np.zeros(shape=(12,36,480,400))
uu = output_data.createVariable('u','float32',\
                               ('ocean_time','s_rho','eta_rho','xi_rho'),\
                               zlib=True,complevel=9,shuffle=True)
vv = output_data.createVariable('v','float32',\
                               ('ocean_time','s_rho','eta_rho','xi_rho'),\
                               zlib=True,complevel=9,shuffle=True)
for i in range(0,12):

  print 'Month %2d' % (i+1)

  #-----------------------------------------------------------------------------
  # Set index limits for current month
  #-----------------------------------------------------------------------------
  stind = (sum(day_count[0:i+1]))*8
  enind = (sum(day_count[0:i+2]))*8+1

  #-----------------------------------------------------------------------------
  # Load and average data for each season
  #-----------------------------------------------------------------------------
  print 'Loading u'
  u =np.zeros(shape=(36,482,401))
  u   =np.nanmean(input_data.variables['u'][stind:enind,:,:,:],axis=0)

  print 'Loading v'
  v =np.zeros(shape=(36,481,402))
  v   =np.nanmean(input_data.variables['v'][stind:enind,:,:,:],axis=0)

  #-----------------------------------------------------------------------------
  # Average velocities to rho-points
  #-----------------------------------------------------------------------------
  uu[i,:,:,:]=(u[:,1:481,1:401]+u[:,1:481,0:400])/2.0
  vv[i,:,:,:]=(v[:,1:481,1:401]+v[:,0:480,1:401])/2.0


input_data.close()
output_data.close()