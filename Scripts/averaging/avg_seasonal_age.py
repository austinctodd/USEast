#===============================================================================
#
# SCRIPT :  avg_seasonal_age.py
#
# PURPOSE : Ingest input data from US East water age model and average the u & v
#           velocities and mean water age for 3 different tracers over 4 seasons.
#           Each variable will be averaged for JFM, AMJ, JAS, and OND.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 16 JANUARY 2015.
#
#===============================================================================

#-------------------------------------------------------------------------------
# Define analysis variables to be updated
#-------------------------------------------------------------------------------
analysis_variables = ('u','v','mean_age_01','mean_age_02','mean_age_03')

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print ""
print "RUNNING: Script avg_seasonal_age.py"
print "----------------------------------------------------------------------"
print ""

##########################################################################
# USER: DO NOT MAKE CHANGES BELOW (if you do, you're on your own!)
##########################################################################

#-------------------------------------------------------------------------------
# Define all required libraries, routines, and modules
#-------------------------------------------------------------------------------
import numpy
import netCDF4
import matplotlib.pyplot as plt
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
input_file  = sys.argv[1]
output_file = sys.argv[2]

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print 'Ingesting data from file %s ...' % (input_file)

#-------------------------------------------------------------------------------
# Open files for reading/writing
#-------------------------------------------------------------------------------
ageanal = Dataset(input_file ,'r')
outanal = Dataset(output_file,'w',format='NETCDF4')

#-------------------------------------------------------------------------------
# Define array dimension attributes
#-------------------------------------------------------------------------------
ageanal_dims = ageanal.dimensions.keys()
outanal_dims = outanal.dimensions.keys()

#-------------------------------------------------------------------------------
# Define array variable attributes
#-------------------------------------------------------------------------------
ageanal_vars = ageanal.variables.keys()
outanal_vars = outanal.variables.keys()

#-------------------------------------------------------------------------------
# Create dimensions in output file
#-------------------------------------------------------------------------------
outanal.createDimension('ocean_time',4)
outanal.createDimension('xi_rho' ,len(ageanal.dimensions['xi_rho' ]))
outanal.createDimension('xi_u'   ,len(ageanal.dimensions['xi_u'   ]))
outanal.createDimension('xi_v'   ,len(ageanal.dimensions['xi_v'   ]))
outanal.createDimension('eta_rho',len(ageanal.dimensions['eta_rho']))
outanal.createDimension('eta_u'  ,len(ageanal.dimensions['eta_u'  ]))
outanal.createDimension('eta_v'  ,len(ageanal.dimensions['eta_v'  ]))
outanal.createDimension('s_rho'  ,len(ageanal.dimensions['s_rho'  ]))

#-------------------------------------------------------------------------------
# Create variables in output file
#-------------------------------------------------------------------------------
lon_rho  = outanal.createVariable('lon_rho',ageanal.variables['lon_rho'].dtype,
                                 (ageanal.variables['lon_rho'].dimensions[0],
                                  ageanal.variables['lon_rho'].dimensions[1]))
lat_rho  = outanal.createVariable('lat_rho',ageanal.variables['lat_rho'].dtype,
                                 (ageanal.variables['lat_rho'].dimensions[0],
                                  ageanal.variables['lat_rho'].dimensions[1]))
#lon_rho[:] = ageanal.variables['salt'][0,34,:,:]
#lat_rho[:] = ageanal.variables['temp'][0,34,:,:]

#var = ageanal.variables['mean_age_01'][0,:,:,:]
#var2 = var.mean(0)
#print var2.shape

#-------------------------------------------------------------------------------
# Loop through variables
#-------------------------------------------------------------------------------
day_count=[31,28,31,30,31,30,31,31,30,31,30,31]

for var in analysis_variables:
  print 'Averaging Variable : %s' % var
  
  #-----------------------------------------------------------------------------
  # Create variable in output file
  #-----------------------------------------------------------------------------
  outvar  = outanal.createVariable(var,ageanal.variables[var].dtype,
                                  (ageanal.variables[var].dimensions[0],
                                   ageanal.variables[var].dimensions[1],
                                   ageanal.variables[var].dimensions[2],
                                   ageanal.variables[var].dimensions[3]))

  #-----------------------------------------------------------------------------
  #  Now average over each season individually
  #-----------------------------------------------------------------------------
  season = 1;
  while (season < 2):
    sys.stdout.write('Season %1i ' % season)

    depth=0
    while depth < 35:
      stind = (sum(day_count[0:(season-1)*3]))*8
      enind = (sum(day_count[0:(season  )*3]))*8-1

      invar = ageanal.variables[var][stind:enind,depth,:,:]
      numpy.mean(invar,axis=0,out=outvar[0,depth,:,:],keepdims=False)
      depth=depth+1

    season=season+1

  print ' '
  print ' '

#print variable_memoryorder

outanal.close()
ageanal.close()


