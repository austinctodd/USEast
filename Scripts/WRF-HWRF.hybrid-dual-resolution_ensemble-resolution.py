#--------------------------------------------------------------------------
#
# SCRIPT:  WRF-HWRF.hybrid-dual-resolution_ensemble-resolution.py
#
# PURPOSE: Ingest a Weather Research and Forecasting (WRF) model Non-
#          hydrostatic Mesoscale Model (NMM) or Hurricane WRF (HWRF)
#          analysis produced by the Gridpoint Statistical Interpolation
#          (GSI) Hybrid Data Assimilation system which has incorporated an
#          ensemble at one resolution and a background at another (i.e.,
#          dual-resolution); an external NETwork Common Data Format
#          (netcdf) file will be created containing the respective
#          prognostic variables from the NMM/HWRF GSI Hybrid analysis at
#          the resolution of the ensemble; this step should proceed the
#          recentering of the ensemble members with respect to the GSI
#          produced hybrid analysis
#
# METHOD:  (1) Define environment variables for data sets to process
#          (2) Run in ~/scripts directory
#
# HISTORY: 06/01/12  Original version. Henry R. Winterbottom
#
#-------------------------------------------------------------------------

# Define analysis variables to be updated

analysis_variables = ('PD','U','V','T','Q','CWM')

# Print message to user

print ""
print "RUNNING: Script WRF-HWRF.hybrid-dual-resolution_ensemble-resolution.py"
print "----------------------------------------------------------------------"
print ""

##########################################################################
# USER: DO NOT MAKE CHANGES BELOW (if you do, you're on your own!)
##########################################################################

# Define all required libraries, routines, and modules

import numpy
import netCDF4
import matplotlib.pyplot as plt
import sys
import mpl_toolkits.basemap
import os

# Define all required sub-libraries, sub-routines, and sub-modules

from netCDF4 import Dataset
from mpl_toolkits.basemap import interp
from numpy import append, arange, dtype, linspace, meshgrid, ndarray

#----

# Get command line arguments

hybridanal_file = sys.argv[1]
ensmemanal_file = sys.argv[2]

# Print message to user

print 'Ingesting variables from GSI hybrid analysis file %s and ENKF ensemble analysis file %s ...' % (hybridanal_file,ensmemanal_file)

# Perform necessary tasks on file

hybridanal = Dataset(hybridanal_file,'r')
ensmemanal = Dataset(ensmemanal_file,'r')
interpanal = Dataset('hybridanalysis_downscale.nc','w',format='NETCDF3_64BIT')

# Define array dimension attributes

hybridanal_dims = hybridanal.dimensions.keys()
ensmemanal_dims = ensmemanal.dimensions.keys()

# Define variable attributes

hybridanal_vars = hybridanal.variables.keys()
ensmemanal_vars = ensmemanal.variables.keys()

# Define local variables

variableattributes = ('FieldType','MemoryOrder','description','units','stagger','coordinates')

#----

# Print message to user

print 'Defining variable dimensions for interpolated hybrid analysis variables in file hybridanalysis_downscale.nc ...'

# Define local variables

interpanal_time           = len(ensmemanal.dimensions['Time'])
interpanal_datestrlen     = len(ensmemanal.dimensions['DateStrLen'])
interpanal_westeast       = len(ensmemanal.dimensions['west_east'])
interpanal_southnorth     = len(ensmemanal.dimensions['south_north'])
interpanal_bottomtop      = len(ensmemanal.dimensions['bottom_top'])
interpanal_bottomtopstag  = len(ensmemanal.dimensions['bottom_top_stag'])
interpanal_soillayersstag = len(ensmemanal.dimensions['soil_layers_stag'])
interpanal_ensemblestag   = len(ensmemanal.dimensions['ensemble_stag'])
interpanal_DIM0006        = len(ensmemanal.dimensions['DIM0006'])
interpanal_DIM0009        = len(ensmemanal.dimensions['DIM0009'])

# Print message to user

print 'Creating array dimension netcdf block within external file hybridanalysis_downscale.nc ...'

# Perform necessary tasks on file

interpanal.createDimension('Time',interpanal_time)
interpanal.createDimension('DateStrLen',interpanal_datestrlen)
interpanal.createDimension('west_east',interpanal_westeast)
interpanal.createDimension('south_north',interpanal_southnorth)
interpanal.createDimension('bottom_top',interpanal_bottomtop)
interpanal.createDimension('bottom_top_stag',interpanal_bottomtopstag)
interpanal.createDimension('soil_layers_stag',interpanal_soillayersstag)
interpanal.createDimension('ensemble_stag',interpanal_ensemblestag)
interpanal.createDimension('DIM0006',interpanal_DIM0006)
interpanal.createDimension('DIM0009',interpanal_DIM0009)

#----

# Print message to user

print 'Defining hybrid analysis geographical grid and preparing for downscaling ... '

# Compute local variables

hybrid_variable_xcoord                        = numpy.linspace(0.0,1.0,len(hybridanal.dimensions['west_east']))
hybrid_variable_ycoord                        = numpy.linspace(0.0,1.0,len(hybridanal.dimensions['south_north']))
ensmem_variable_xcoord                        = numpy.linspace(0.0,1.0,len(ensmemanal.dimensions['west_east']))
ensmem_variable_ycoord                        = numpy.linspace(0.0,1.0,len(ensmemanal.dimensions['south_north']))
ensmem_variable_xcoord,ensmem_variable_ycoord = numpy.meshgrid(ensmem_variable_xcoord,ensmem_variable_ycoord)

#----

# Loop through each analysis variable specified by user and downscale
# respective variables to lower-resolution analysis

for analysis_variable in analysis_variables:

    # Define local variables

    hybrid_variable_memoryorder = getattr(hybridanal.variables[analysis_variable],'MemoryOrder')

    # Define local variables and compute output variables accordingly

    if(hybrid_variable_memoryorder == 'XY '):

        # Print message to user

        print 'Ingesting 2-dimensional variable %s from file %s ...' % (analysis_variable,hybridanal_file)

        # Define local variable

        hybrid_variable_array = hybridanal.variables[analysis_variable][0,:,:]

        # Define local variables

        analysis_variable_strname_tdim = ensmemanal.variables[analysis_variable].dimensions[0]
        analysis_variable_strname_ydim = ensmemanal.variables[analysis_variable].dimensions[1]
        analysis_variable_strname_xdim = ensmemanal.variables[analysis_variable].dimensions[2]

        # Initialize local variable

        analysis_variable_array = ensmemanal.variables[analysis_variable][:,:,:]

        # Print message to user

        print 'Interpolating variable %s from horizontal resolution within %s to horizontal resolution of %s ...' % (analysis_variable,hybridanal_file,ensmemanal_file)

        # Compute local variable

        analysis_variable_array = interp(hybrid_variable_array,hybrid_variable_xcoord,hybrid_variable_ycoord,ensmem_variable_xcoord,ensmem_variable_ycoord,order=3)

        # Define local variable

        analysis_variable_2d = analysis_variable_array

        # Perform necessary tasks on file

        analysis_variable_2d    = interpanal.createVariable(analysis_variable,ensmemanal.variables[analysis_variable].dtype,(analysis_variable_strname_tdim,analysis_variable_strname_ydim,analysis_variable_strname_xdim))
        analysis_variable_2d[:] = analysis_variable_array

        # Print message to user

        print 'Defining variable attributes for variable %s ...' % (analysis_variable)

        # Loop through each variable attribute, define local
        # variables, and perform necessary tasks on file

        for variableattribute in variableattributes:

            # Define local variables

            attribute_str = '%s' % (variableattribute)

            # Inquire as to whether attribute exists for variable;
            # if so, obtain value and write to external netcdf
            # formatted file

            if(hasattr(ensmemanal.variables[analysis_variable],attribute_str)):

                # Define local variable
                
                attribute_val = getattr(ensmemanal.variables[analysis_variable],attribute_str)

                # Perform necessary tasks on file
            
                setattr(analysis_variable_2d,variableattribute,attribute_val)
                                 
    # Define local variables accordingly and compute output variables
    # accordingly

    if(hybrid_variable_memoryorder == 'XYZ'):

        # Print message to user

        print 'Ingesting 3-dimensional variable %s from file %s ...' % (analysis_variable,hybridanal_file)

        # Define local variables

        analysis_variable_strname_tdim = ensmemanal.variables[analysis_variable].dimensions[0]
        analysis_variable_strname_zdim = ensmemanal.variables[analysis_variable].dimensions[1]
        analysis_variable_strname_ydim = ensmemanal.variables[analysis_variable].dimensions[2]
        analysis_variable_strname_xdim = ensmemanal.variables[analysis_variable].dimensions[3]

        # Define local variable accordingly

        if(analysis_variable_strname_zdim == 'bottom_top'):

            # Define local variable

            numlevels = interpanal_bottomtop

        # Define local variable accordingly

        if(analysis_variable_strname_zdim == 'bottom_top_stag'):

            # Define local variable

            numlevels = interpanal_bottomtopstag

        # Print message to user

        print 'Interpolating variable %s from horizontal resolution within %s to horizontal resolution of %s ...' % (analysis_variable,hybridanal_file,ensmemanal_file)

        # Initialize local variable

        local_analysis_variable_array = ensmemanal.variables[analysis_variable][0,:,:,:]

        # Initialize local variable

        level = 0

        # Loop through vertical coordinate and ingest, compute, and
        # update variables from/within external file

        while(level < numlevels):

            # Print message to user

            print 'Interpolating level %i for variable %s ...' % (level,analysis_variable)

            # Define local variables

            local_hybrid_variable_array = hybridanal.variables[analysis_variable][0,level,:,:]

            # Compute local variable

            local_analysis_variable_array = interp(local_hybrid_variable_array,hybrid_variable_xcoord,hybrid_variable_ycoord,ensmem_variable_xcoord,ensmem_variable_ycoord,order=3)

            # Define local variable accordingly

            if(level == 0):

                # Define local variable

                local_analysis_variable_3d = local_analysis_variable_array

            else:

                # Define local variable

                local_analysis_variable_3d = append(local_analysis_variable_3d,local_analysis_variable_array)

            # Update counting variable

            level += 1

        # Perform necessary tasks on file

        analysis_variable_3d    = local_analysis_variable_3d
        analysis_variable_3d    = interpanal.createVariable(analysis_variable,ensmemanal.variables[analysis_variable].dtype,(analysis_variable_strname_tdim,analysis_variable_strname_zdim,analysis_variable_strname_ydim,analysis_variable_strname_xdim))
        analysis_variable_3d[:] = local_analysis_variable_3d

        # Print message to user

        print 'Defining variable attributes for variable %s ...' % (analysis_variable)

        # Loop through each variable attribute, define local
        # variables, and perform necessary tasks on file

        for variableattribute in variableattributes:

            # Define local variables

            attribute_str = '%s' % (variableattribute)

            # Inquire as to whether attribute exists for variable;
            # if so, obtain value and write to external netcdf
            # formatted file

            if(hasattr(ensmemanal.variables[analysis_variable],attribute_str)):

                # Define local variable
                
                attribute_val = getattr(ensmemanal.variables[analysis_variable],attribute_str)

                # Perform necessary tasks on file
            
                setattr(analysis_variable_3d,variableattribute,attribute_val)

#----

# Print message to user

print 'Closing all external files ingested and created by routine ...'

# Perform necessary tasks on file

hybridanal.close()
ensmemanal.close()
interpanal.close()
