#===============================================================================
#
# SCRIPT :  add_fwd_vars.py
#
# PURPOSE : Ingest data from US East water age forward model and average over
#           the five years of forward model state.  Several variables already
#           have been averaged over the 5 year period for use in the LTRANS
#           simulations.  Several new netCDF dimensions are also added to the
#           existing file.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 11 February 2015.
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
print " RUNNING: Script add_fwd_vars.py"
print "----------------------------------------------------------------------"
print ""

################################################################################
#                                                                              #
#                            User Defined Variables                            #
#                                                                              #
infile1 = '/gpfs_common/he_share/actodd/USeast-age/output/clim/fwdrun/part1/'+\
          'useast_his.nc'
infile2 = '/gpfs_common/he_share/actodd/USeast-age/output/clim/fwdrun/part2/'+\
          'useast_his.nc'
outfile = '/gpfs_common/he_share/actodd/USeast-age/output/clim/fwdrun/'+\
           'daily_avg.nc'
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Read in additional grid variables from forward file
#-------------------------------------------------------------------------------
indata1=Dataset(inpfile1,mode='r')
indata2=Dataset(inpfile2,mode='r')

#-------------------------------------------------------------------------------
# Open output file and create dimensions and variables
#-------------------------------------------------------------------------------
print " "
print "Opening output file: %s" % output_file
print " "
outdata=Dataset(outfile,mode='w')

#--------------------------------------------------------------------------
# Declare 2d variables to be averaged
#--------------------------------------------------------------------------
analysis_variables2d = ('DU_avg1','DU_avg2','DV_avg1','DV_avg2','shflux',\
                        'ssflux','swrad','sustr','sustr','svstr')

for analysis_variable in analysis_variables2d:
  print ""
  print "---------------------------------------------------------"
  print " Beginning Averging for variable: %s" % analysis_variable
  print "---------------------------------------------------------"
  print ""

  #--------------------------------------------------------------------------
  # Define variable in output dataset
  #--------------------------------------------------------------------------
  var=outdata.createVariable(analysis_variable,indata1.variables[analysis_variable].dtype,\
                             (indata1.dimensions[analysis_variable][0],\
                              indata1.dimensions[analysis_variable][1],\
                              indata1.dimensions[analysis_variable][2]))

  #--------------------------------------------------------------------------
  # Loop through the rest of the times
  #--------------------------------------------------------------------------
  print ""
  print "---------------------------------------------------"
  print "           Beginning Main Loop Iteration           "
  print "---------------------------------------------------"
  print ""

  tdim=0
  while tdim < 2923:

    print 'Time : %04i / 2923' % (tdim+1)
    #------------------------------------------------------------------------
    # Load in the SSH, Temp, and Salt to calculate layer & thermocline depths
    #------------------------------------------------------------------------
    inds=[0,2922,0,2922,5844,8766]+tdim
    var[tdim,:,:]=(indata1.variables[analysis_variable][inds[0],:,:]+\
                   indata1.variables[analysis_variable][inds[1],:,:]+\
                   indata2.variables[analysis_variable][inds[2],:,:]+\
                   indata2.variables[analysis_variable][inds[3],:,:]+\
                   indata2.variables[analysis_variable][inds[4],:,:]+\
                   indata2.variables[analysis_variable][inds[5],:,:])/6.0

     tdim=tdim+1

#--------------------------------------------------------------------------
# Now average 3D variables
#--------------------------------------------------------------------------
analysis_variables3d = ('AKv','AKt')
for analysis_variable in analysis_variables2d:
  print ""
  print "---------------------------------------------------------"
  print " Beginning Averging for variable: %s" % analysis_variable
  print "---------------------------------------------------------"
  print ""
  
  #--------------------------------------------------------------------------
  # Define variable in output dataset
  #--------------------------------------------------------------------------
  var=outdata.createVariable(analysis_variable,indata1.variables[analysis_variable].dtype,\
                             (indata1.dimensions[analysis_variable][0],\
                              indata1.dimensions[analysis_variable][1],\
                              indata1.dimensions[analysis_variable][2])\
                              indata1.dimensions[analysis_variable][3]))
                              
  #--------------------------------------------------------------------------
  # Loop through the rest of the times
  #--------------------------------------------------------------------------
  print ""
  print "---------------------------------------------------"
  print "           Beginning Main Loop Iteration           "
  print "---------------------------------------------------"
  print ""
                              
  tdim=0
  while tdim < 2923:
                                
    print 'Time : %04i / 2923' % (tdim+1)

    #------------------------------------------------------------------------
    # Load in the SSH, Temp, and Salt to calculate layer & thermocline depths
    #------------------------------------------------------------------------
    inds=[0,2922,0,2922,5844,8766]+tdim
    var[tdim,:,:,:]=(indata1.variables[analysis_variable][inds[0],:,:,:]+\
                     indata1.variables[analysis_variable][inds[1],:,:,:]+\
                     indata2.variables[analysis_variable][inds[2],:,:,:]+\
                     indata2.variables[analysis_variable][inds[3],:,:,:]+\
                     indata2.variables[analysis_variable][inds[4],:,:,:]+\
                     indata2.variables[analysis_variable][inds[5],:,:,:])/6.0
                                      
    tdim=tdim+1

indata1.close()
indata2.close()
outdata.close()




