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
analysis_variables = ('var')

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print ""
print "----------------------------------------------------------------------"
print " RUNNING: Script test.py"
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
import mpl_toolkits.basemap
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
inanal = Dataset(input_file ,'r')
outanal = Dataset(output_file,'w',format='NETCDF4')

#-------------------------------------------------------------------------------
# Define array dimension attributes
#-------------------------------------------------------------------------------
inanal_dims = inanal.dimensions.keys()
outanal_dims = outanal.dimensions.keys()

#-------------------------------------------------------------------------------
# Define array variable attributes
#-------------------------------------------------------------------------------
inanal_vars = inanal.variables.keys()
outanal_vars = outanal.variables.keys()

#-------------------------------------------------------------------------------
# Create dimensions in output file
#-------------------------------------------------------------------------------
outanal.createDimension('x',1)
outanal.createDimension('y',len(inanal.dimensions['y']))

#-------------------------------------------------------------------------------
# Create variables in output file
#-------------------------------------------------------------------------------
var  = numpy.array(inanal.variables['var'])
var2 = var.mean(0)

sys.stdout.write('test')
print var2.shape

outvar = var2
outvar  = outanal.createVariable('var',inanal.variables['var'].dtype,
                                 (inanal.variables['var'].dimensions[0],
                                  inanal.variables['var'].dimensions[1]))
outvar[:] = var2

#var = numpy.array(inanal.variables['mean_age_01'][0,:,:,:]
#var2 = var.mean(0)
#print var2.shape

outanal.close()
inanal.close()


import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset

infile='/Volumes/Black_box/Data/USeast-rtime/output/useast_rtime_0010.nc'

fh=Dataset(infile,mode='r')
lon =fh.variables['lon_rho'  ][:]
lat =fh.variables['lat_rho'  ][:]
age =fh.variables['age_01'   ][100,35,:,:]
mask=fh.variables['scope_rho'][:]
fh.close()

lon_0=lon.mean()
lat_0=lat.mean()
m=Basemap(projection='merc',llcrnrlat=8,urcrnrlat=45,llcrnrlon=-98,urcrnrlon=-60,lat_ts=20,resolution='m')

# Plot Data
xi, yi = m(lon,lat)
cs = m.pcolor(xi,yi,np.squeeze(age))

# Add Grid Lines
m.drawparallels(np.arange(0.,90.,10.))
m.drawmeridians(np.arange(-100.,-40.,10.))

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label('residence time')

# Add Title
plt.title('Age tracer 1')

plt.show()

