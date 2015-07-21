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
import sys
import os

#-------------------------------------------------------------------------------
# Define all required sub-libraries, sub-routines, and sub-modules
#-------------------------------------------------------------------------------
from netCDF4 import Dataset

#-------------------------------------------------------------------------------
# Print message to screen
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script avg_transit_times.py"
print "----------------------------------------------------------------------"
print ""

################################################################################
#                                                                              #
#                            User Defined Variables                            #
#                                                                              #
age_file    = '/gpfs_share/actodd/USeast-age/output/clim/daily_avg.nc'
ATL_dir     = '/gpfs_share/actodd/USeast-rtime/input/ATL/'
GOM_dir     = '/gpfs_share/actodd/USeast-rtime/input/GOM_shelf/'
grid_fileA  = '/gpfs_share/actodd/USeast-rtime/input/In/grid_ATLscope.nc'
grid_fileG  = '/gpfs_share/actodd/USeast-rtime/input/In/grid_GOM_shelf_scope.nc'
output_file = '/gpfs_common/he_share/actodd/USeast-age/output/clim/'+\
              'avg_transit_times.nc'
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Open grid file and read in variables
#-------------------------------------------------------------------------------
grid_data=Dataset(grid_fileA,mode='r')
hmask  = grid_data.variables['mask_rho' ][:]
Ascope = grid_data.variables['scope_rho'][:]
mask   = grid_data.variables['mask_rho' ][:,:]
lon    = grid_data.variables['lon_rho'  ][:,:]
lat    = grid_data.variables['lat_rho'  ][:,:]
h      = grid_data.variables['h'        ][:,:]
pm     = grid_data.variables['pm'       ][:]
pn     = grid_data.variables['pn'       ][:]
grid_data.close()

# Get GOM Shelf gird
grid_data=Dataset(grid_fileG,mode='r')
Gscope = grid_data.variables['scope_rho'][:]
grid_data.close()

hc=200.
N=36

#-------------------------------------------------------------------------------
# Read in additional grid variables from forward file
#-------------------------------------------------------------------------------
age_data=Dataset(age_file,mode='r')
s_r =age_data.variables['s_rho'][:]
s_w =age_data.variables['s_w'  ][:]
Cs_r=age_data.variables['Cs_r' ][:]
Cs_w=age_data.variables['Cs_w' ][:]

L=len(age_data.dimensions['xi_rho' ])
M=len(age_data.dimensions['eta_rho'])
N=len(age_data.dimensions['s_rho'  ])

#-------------------------------------------------------------------------------
# Open output file and create dimensions and variables
#-------------------------------------------------------------------------------
print " "
print "Creating output file: %s" % output_file
print " "

output_data=Dataset(output_file,mode='w',format='NETCDF4')
output_data.createDimension('ocean_time',12)
output_data.createDimension('s_rho'  ,len(age_data.dimensions['s_rho'  ]))
output_data.createDimension('xi_rho' ,len(age_data.dimensions['xi_rho' ]))
output_data.createDimension('eta_rho',len(age_data.dimensions['eta_rho']))

lon_rho = output_data.createVariable('lon_rho','float32',('eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
lat_rho = output_data.createVariable('lat_rho','float32',('eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)
#--------------------------------------------------------------------------
#  Create standard variables (lon,lat,etc.)
#--------------------------------------------------------------------------
lon_rho[ :] = age_data.variables['lon_rho'   ][:]
lat_rho[ :] = age_data.variables['lat_rho'   ][:]

#--------------------------------------------------------------------------
# Loop through the water age output to calculate monthly means
#--------------------------------------------------------------------------
print "------------------------------"
print "      Beginning Age Loop      "
print "------------------------------"

day_count=[0,31,28,31,30,31,30,31,31,30,31,30,31]
analysis_variables = ('mean_age_01','mean_age_02','mean_age_03')
for analysis_variable in analysis_variables:
    
  print 'Variable: %s' % analysis_variable

  var = np.zeros(shape=(12,N,M,L))
  var = output_data.createVariable(analysis_variable,'float32',\
                                  ('ocean_time','s_rho','eta_rho','xi_rho'),\
                                   zlib=True,complevel=9,shuffle=True)
  for tdim in range(0,12):

    print 'Month : %04i' % (tdim+1)

    #------------------------------------------------------------------------
    # Set index limits for current month
    #------------------------------------------------------------------------
    stind = (sum(day_count[0:tdim+1]))*8
    enind = (sum(day_count[0:tdim+2]))*8
    
    #------------------------------------------------------------------------
    # ATL mean residence times
    #------------------------------------------------------------------------
    tvar =age_data.variables[analysis_variable][stind:enind,:,:,:]
    tvar =np.nanmean(tvar,axis=0)
    var[tdim,:,:,:]=tvar

age_data.close()
output_data.close()




