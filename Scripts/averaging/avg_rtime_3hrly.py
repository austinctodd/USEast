#===============================================================================
#
# SCRIPT :  avg_rtime_3hrly.py
#
# PURPOSE : Ingest input data from US East residence time model and average the 
#           calculated residence time variables over the 5 analysis years for
#           every 3 hours.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 20 JANUARY 2015.
#
#===============================================================================

#-------------------------------------------------------------------------------
# Define analysis variables to be updated
#-------------------------------------------------------------------------------
analysis_variables = ('age_01','age_02')

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print ""
print "----------------------------------------------------------------------"
print "RUNNING: Script avg_rtime_3hrly.py"
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
#import matplotlib.pyplot as plt
import sys
#import mpl_toolkits.basemap
import os

#-------------------------------------------------------------------------------
# Define all required sub-libraries, sub-routines, and sub-modules
#-------------------------------------------------------------------------------
from netCDF4 import Dataset
#from mpl_toolkits.basemap import interp
from numpy import append, arange, dtype, linspace, meshgrid, ndarray

#-------------------------------------------------------------------------------
# Get command line arguments
#-------------------------------------------------------------------------------
input_files = ('useast_rtime_0010.nc','useast_rtime_0009.nc',
               'useast_rtime_0008.nc','useast_rtime_0007.nc',
	       'useast_rtime_0006.nc');
output_file = 'rtime_avg_3hrly.nc'

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print 'Averaging data to file %s ...' % (output_file)
print 'First file: %s' % input_files[0]

#-------------------------------------------------------------------------------
# Open files for reading/writing
#-------------------------------------------------------------------------------
year1anal = Dataset(input_files[0] ,'r')
year2anal = Dataset(input_files[1] ,'r')
year3anal = Dataset(input_files[2] ,'r')
year4anal = Dataset(input_files[3] ,'r')
year5anal = Dataset(input_files[4] ,'r')
outpuanal = Dataset(output_file,'w',format='NETCDF4')

#-------------------------------------------------------------------------------
# Define array dimension attributes
#-------------------------------------------------------------------------------
year1anal_dims = year1anal.dimensions.keys()
year2anal_dims = year2anal.dimensions.keys()
year3anal_dims = year3anal.dimensions.keys()
year4anal_dims = year4anal.dimensions.keys()
year5anal_dims = year5anal.dimensions.keys()
outpuanal_dims = outpuanal.dimensions.keys()

#-------------------------------------------------------------------------------
# Define array variable attributes
#-------------------------------------------------------------------------------
year1anal_vars = year1anal.variables.keys()
year2anal_vars = year2anal.variables.keys()
year3anal_vars = year3anal.variables.keys()
year4anal_vars = year4anal.variables.keys()
year5anal_vars = year5anal.variables.keys()
outpuanal_vars = outpuanal.variables.keys()

#-------------------------------------------------------------------------------
# Create dimensions in output file
#-------------------------------------------------------------------------------
outpuanal.createDimension('ocean_time',len(year1anal.dimensions['ocean_time']))
outpuanal.createDimension('xi_rho'    ,len(year1anal.dimensions['xi_rho'    ]))
outpuanal.createDimension('eta_rho'   ,len(year1anal.dimensions['eta_rho'   ]))
outpuanal.createDimension('s_rho'     ,len(year1anal.dimensions['s_rho'     ]))

#-------------------------------------------------------------------------------
# Create variables in output file
#-------------------------------------------------------------------------------
lon_rho = outpuanal.createVariable('lon_rho',year1anal.variables['lon_rho'].dtype,
                                  (year1anal.variables['lon_rho'].dimensions[0],
                                  year1anal.variables['lon_rho'].dimensions[1]))
lat_rho = outpuanal.createVariable('lat_rho',year1anal.variables['lat_rho'].dtype,
                                 (year1anal.variables['lat_rho'].dimensions[0],
                                  year1anal.variables['lat_rho'].dimensions[1]))
age_01  = outpuanal.createVariable('age_01',year1anal.variables['age_01'].dtype,
                                 (year1anal.variables['age_01'].dimensions[0],
                                  year1anal.variables['age_01'].dimensions[1],
                                  year1anal.variables['age_01'].dimensions[2],
                                  year1anal.variables['age_01'].dimensions[3]))
age_02  = outpuanal.createVariable('age_02',year1anal.variables['age_02'].dtype,
                                 (year1anal.variables['age_02'].dimensions[0],
                                  year1anal.variables['age_02'].dimensions[1],
                                  year1anal.variables['age_02'].dimensions[2],
                                  year1anal.variables['age_02'].dimensions[3]))
lon_rho[:] = year1anal.variables['lon_rho']
lat_rho[:] = year1anal.variables['lat_rho']

#var = year1anal.variables['mean_age_01'][0,:,:,:]
#var2 = var.mean(0)
#print var2.shape

#-------------------------------------------------------------------------------
# Loop through all times
#-------------------------------------------------------------------------------
#day_count=[31,28,31,30,31,30,31,31,30,31,30,31]

tdim = 0 # Time dimension
while tdim < len(year1anal.dimensions['ocean_time']):

  print 'Time %04i' % (tdim)

  srho = 0
  while srho < 35:
    
    sys.stdout.write('%02i.' % srho)
  
    a1 = np.array(year1anal.variables['age_01'][tdim,srho,:,:])
    a2 = np.array(year2anal.variables['age_01'][tdim,srho,:,:])
    a3 = np.array(year3anal.variables['age_01'][tdim,srho,:,:])
    a4 = np.array(year4anal.variables['age_01'][tdim,srho,:,:])
    a5 = np.array(year5anal.variables['age_01'][tdim,srho,:,:])
    age_01[tdim,srho,:,:] = (a1 + a2 + a3 + a4 + a5) / 5.0		  

    a1 = np.array(year1anal.variables['age_02'][tdim,srho,:,:])
    a2 = np.array(year2anal.variables['age_02'][tdim,srho,:,:])
    a3 = np.array(year3anal.variables['age_02'][tdim,srho,:,:])
    a4 = np.array(year4anal.variables['age_02'][tdim,srho,:,:])
    a5 = np.array(year5anal.variables['age_02'][tdim,srho,:,:])
    age_02[tdim,srho,:,:] = (a1 + a2 + a3 + a4 + a5) / 5.0		  

    srho=srho+1

  tdim=tdim+1

#for var in analysis_variables:
#  print 'Averaging Variable : %s' % var
  
  #-----------------------------------------------------------------------------
  # Create variable in output file
  #-----------------------------------------------------------------------------
#  outvar  = outanal.createVariable(var,year1anal.variables[var].dtype,
#                                  (year1anal.variables[var].dimensions[0],
#                                   year1anal.variables[var].dimensions[1],
#                                   year1anal.variables[var].dimensions[2],
#                                   year1anal.variables[var].dimensions[3]))

  #-----------------------------------------------------------------------------
  #  Now average over each season individually
  #-----------------------------------------------------------------------------
#  season = 1;
#  while (season < 2):
#    sys.stdout.write('Season %1i ' % season)
#
#    depth=0
#    while depth < 35:
#      stind = (sum(day_count[0:(season-1)*3]))*8
#      enind = (sum(day_count[0:(season  )*3]))*8-1
#
#      invar = year1anal.variables[var][stind:enind,depth,:,:]
#      numpy.mean(invar,axis=0,out=outvar[0,depth,:,:],keepdims=False)
#      depth=depth+1
#
#    season=season+1
#
#  print ' '
#  print ' '
#
#print variable_memoryorder

outanal.close()
year1anal.close()
year2anal.close()
year3anal.close()
year4anal.close()
year5anal.close()


