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

#-------------------------------------------------------------------------------
# Print message to screen
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script transit_time_above_pcline.py"
print "----------------------------------------------------------------------"
print ""

################################################################################
#                                                                              #
#                            User Defined Variables                            #
#                                                                              #
input_file  = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
              'avg_3hrly.nc'
Grtime_file = '/Volumes/Black_box/Data/USeast-rtime/output/GOM/'+\
              'rtime_GOM.nc'
GOM_file    = '/Volumes/Black_box/Data/USeast/Data/grd/grid_GOMscope.nc'
grid_file   = '/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc'
output_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
              'pcline_avgs.nc'
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Open grid file and read in variables
#-------------------------------------------------------------------------------
print 'Ingesting data from file %s ...' % (grid_file)
grid_data=Dataset(grid_file,mode='r')
lat =grid_data.variables['lat_rho'  ][:,:]
h   =grid_data.variables['h'        ][:,:]
grid_data.close()

# Get GOM gird
print 'Ingesting data from file %s ...' % (GOM_file)
grid_data=Dataset(GOM_file,mode='r')
Gscope = grid_data.variables['scope_rho'][:]
grid_data.close()

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

hc=200.

#-------------------------------------------------------------------------------
# Open residence time files
#-------------------------------------------------------------------------------
Gtime_data =Dataset(Grtime_file,mode='r')

#--------------------------------------------------------------------------
#  Create local storage variables
#--------------------------------------------------------------------------
rtime01 = np.zeros(shape=(len(input_data.variables['ocean_time']),M,L))

#--------------------------------------------------------------------------
# Loop through the water age output to calculate thermocline depth
#--------------------------------------------------------------------------
print "-----------------------------"
print "Beginning Main Loop Iteration"
print "-----------------------------"
for tdim in range(0,len(input_data.variables['ocean_time'])):

  print 'Time : %04i' % tdim

  for i in range(0,M):
    for j in range(0,L):
      
      # Cell is over the water
      if Gscope[i,j]>0:

        #-----------------------------------------------------------------
        # Load in temp,salt, age, residence time
        #----------------------------------------------------------------
        zeta  =input_data.variables['zeta'  ][tdim,  i,j]
        temp  =input_data.variables['temp'  ][tdim,:,i,j]
        salt  =input_data.variables['salt'  ][tdim,:,i,j]
        Grtime=Gtime_data.variables['age_02'][tdim,:,i,j]

        #------------------------------------------------------------------------
        # Calculate depths of each rho layer
        #------------------------------------------------------------------------
        z =np.zeros(shape=(N,))
        zw=np.zeros(shape=(N+1,))
        for k in range(0,len(s_r)):
          z0  =(hc*s_r[k]+Cs_r[k]*h[i,j])/(hc+h[i,j]);
          z[k]=zeta+(zeta+h[i,j])*z0;

          z0   =(hc*s_w[k]+Cs_w[k]*h[i,j])/(hc+h[i,j]);
          zw[k]=zeta+(zeta+h[i,j])*z0;

        # Add last depth for zw
        z0   =(hc*s_w[N]+Cs_w[N]*h[i,j])/(hc+h[i,j]);
        zw[N]=zeta+(zeta+h[i,j])*z0;

        Hz = abs(zw[1:N+1]-zw[0:N])
  
        #------------------------------------------------------------------------
        # Calculate various seawater properties
        #------------------------------------------------------------------------
        pres=sw.pres(z*-1,lat[i,j])

        #-----------------------------------------------------------------
        # Calculate depth of the pycnocline
        #----------------------------------------------------------------
        dens=sw.dens(salt,temp,pres)-1000.
        drdZ = abs(dens[1:N]-dens[0:N-1])/abs(z[1:N]-z[0:N-1])
        drdZ_ind=np.argmax(drdZ)
        
        #-----------------------------------------------------------------
        # Average residence time above the pycnocline
        #----------------------------------------------------------------
        rtime01[tdim,i,j]=sum(   Hz[drdZ_ind::]*Grtime[drdZ_ind::])/\
                              abs(z[drdZ_ind])
      else:
        rtime01[tdim,i,j]=np.nan
    
#-------------------------------------------------------------------------------
# Open output file and create dimensions and variables
#-------------------------------------------------------------------------------
print " "
print "Creating output file: %s" % output_file
print " "

output_data=Dataset(output_file,mode='a')
rtime_1 = output_data.createVariable('GOM_rtime2','float32',\
                                     ('ocean_time','eta_rho','xi_rho'),\
                                     zlib=True,complevel=9,shuffle=True)

#--------------------------------------------------------------------------
#  Write output data to file
#--------------------------------------------------------------------------
rtime_1[:]  = rtime01[:]

# Close Files
print 'Closing files'
input_data.close()
Gtime_data.close()
output_data.close()




