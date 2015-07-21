#===============================================================================
#
# SCRIPT :  plot_depav_vels_shelf.py
#
# PURPOSE : Ingest input data from US East water age model and plot the monthly
#           mean velocities on a map.  Try masking out non-shelf values
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 1 JULY 2015 for personal use.
#
#===============================================================================

input_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
             'monthly_avg_vels.nc'
ATL_file   = '/Volumes/Black_box/Data/USeast/Data/grd/grid_ATLscope.nc'
fwd_file   = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/'+\
              'avg_3hrly.nc'
plot_dir   = '/Volumes/Black_box/Data/PLOTS/USeast-age/vels/ATL/depav/'+\
             'shelf/transparent/'

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script plot_mean_vels.py"
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
#import set_depth
#import PIL

#-------------------------------------------------------------------------------
# Define all required sub-libraries, sub-routines, and sub-modules
#-------------------------------------------------------------------------------
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
#from numpy import append, arange, dtype, linspace, meshgrid, ndarray

#-------------------------------------------------------------------------------
# Open grid file and read in variables
#-------------------------------------------------------------------------------
print 'Ingesting data from file %s ...' % (ATL_file)
grid_data=Dataset(ATL_file,mode='r')
hmask = grid_data.variables['mask_rho' ][:,1:401-24]
scope = grid_data.variables['scope_rho'][:,1:401-24]
h     = grid_data.variables['h'        ][1:481,1:401-24]
lon   = grid_data.variables['lon_rho'  ][:,1:401-24]
lat   = grid_data.variables['lat_rho'  ][:,1:401-24]
grid_data.close()

#-------------------------------------------------------------------------------
# Open transit time and fwd file
#-------------------------------------------------------------------------------
input_data=Dataset(input_file,mode='r')
lonu = lon[1:481,:]
latu = lat[1:481,:]

lon2=lonu-0.5*(lonu[1,2]-lonu[1,1])
lat2=latu-0.5*(latu[2,1]-latu[1,1])
print 'min lon %d max lon %d' % (lon2.min(),lon2.max())

#-------------------------------------------------------------------------------
# Read in additional grid variables from forward file
#-------------------------------------------------------------------------------
fwd_data=Dataset(fwd_file,mode='r')
s_r =fwd_data.variables['s_rho'][:]
s_w =fwd_data.variables['s_w'  ][:]
Cs_r=fwd_data.variables['Cs_r' ][:]
Cs_w=fwd_data.variables['Cs_w' ][:]

L=len(fwd_data.dimensions['xi_rho' ])-24
M=len(fwd_data.dimensions['eta_rho'])
N=len(fwd_data.dimensions['s_rho'  ])

hc=200.

#-------------------------------------------------------------------------------
# Mask the age variables to the shelf
#-------------------------------------------------------------------------------
scope2=scope+hmask
scope=1-scope
hmask=1-hmask
cm.gist_ncar_r.set_bad('white',alpha=0)

#-------------------------------------------------------------------------------
# Open age file
#-------------------------------------------------------------------------------
day_count=[0,31,28,31,30,31,30,31,31,30,31,30,31]
fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor='w')
cm.spectral_r.set_bad('grey',alpha=0)

for i in range(0,1):

  print 'Month %2d' % (i+1)

  #-----------------------------------------------------------------------------
  # Load and average velocities for each season
  #-----------------------------------------------------------------------------
  print 'Loading u'
  #u   =np.nanmean(input_data.variables['u'][i:i+3,0,:,:],axis=0)
  #v   =np.nanmean(input_data.variables['v'][i:i+3,0,:,:],axis=0)
  u   =input_data.variables['u'][i,:,:,0:401-25]
  v   =input_data.variables['v'][i,:,:,0:401-25]
  umag=np.sqrt(np.square(u) + np.square(v))

  #------------------------------------------------------------------------
  # Load in temp,salt, age, residence time
  #------------------------------------------------------------------------
  stind=sum(day_count[:i+1])*8
  enind=sum(day_count[:i+2])*8
  zeta  =np.nanmean(fwd_data.variables['zeta'  ][stind:enind,1:481,1:401-24],axis=0)
  
  #------------------------------------------------------------------------
  # Calculate depths of each rho layer
  #------------------------------------------------------------------------
  z =np.zeros(shape=(N  ,M-2,L-2))
  zw=np.zeros(shape=(N+1,M-2,L-2))
  for k in range(0,len(s_r)):
    z0  =(hc*s_r[k]+Cs_r[k]*h)/(hc+h);
    z[k,:,:]=zeta+(zeta+h)*z0;
    
    z0   =(hc*s_w[k]+Cs_w[k]*h)/(hc+h);
    zw[k,:,:]=zeta+(zeta+h)*z0;

  # Add last depth for zw
  z0   =(hc*s_w[N]+Cs_w[N]*h)/(hc+h);
  zw[N,:,:]=zeta+(zeta+h)*z0;
  
  Hz = abs(zw[1:N+1,:,:]-zw[0:N,:,:])

  #------------------------------------------------------------------------
  # Multiply velocities by layer depths
  #------------------------------------------------------------------------
  u2    = np.nansum(Hz*u   ,axis=0)/(h+zeta)
  v2    = np.nansum(Hz*v   ,axis=0)/(h+zeta)
  umag2 = np.nansum(Hz*umag,axis=0)/(h+zeta)

  u2   =np.ma.array(u2   ,mask=scope[1:481,:])
  v2   =np.ma.array(v2   ,mask=scope[1:481,:])
  umag2=np.ma.array(umag2,mask=scope[1:481,:])
 
  #-----------------------------------------------------------------------------
  # Plot Data using Basemap mercator projection
  #-----------------------------------------------------------------------------
  print 'Plotting data.'
  plt.clf()

  # Plot Data
  if i==18:
    m=Basemap(projection='merc',llcrnrlat=26.5,urcrnrlat=lat[481,401-26]+8,\
              llcrnrlon=-90.75,urcrnrlon=lon[481,401-26]+5,lat_ts=20,resolution='f')

    # Draw a think border around the whole map
    m.drawmapboundary(linewidth=3)

    # Plot Data
    xi, yi = m(lon2,lat2)
    cs = m.pcolormesh(xi,yi,umag2,cmap=cm.gist_ncar_r,vmin=0.005, vmax = 0.15)
    xi, yi = m(lon,lat)
    cM = m.contour(   xi,yi,scope2,[1,1],colors='k',linewidth=2.0)
  
    # transform vectors to projection grid.
    uproj,vproj,xx,yy = m.transform_vector(u2[1::2,1::2]*7.5,v2[1::2,1::2]*7.5,\
                                           lonu[0,1::2].squeeze(),latu[1::2,0].squeeze(),\
                                           60,60,returnxy=True,masked=True)
    # now plot.
    Q  = m.quiver(xx,yy,uproj,vproj,scale=15)

    # Add Coastlines, States, and Country Boundaries
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.fillcontinents(color='grey')
  
    # Add Colorbar
    clabels=np.array([0.005,0.025,0.05,0.075,0.1,0.12,0.15])#,2.75,3])
    cbar = m.colorbar(cs, location='right', pad="2%",ticks=clabels,extend='max')
    cbar.ax.set_yticklabels([0.5,2.5,5,7.5,10,12,15])#,2.75,3])
    cbar.set_label('Current Speed (cm/s)')

  else:
    fig=plt.figure(num=None, figsize=(10, 10), dpi=300, facecolor=None)
    plt.clf()
    fig.patch.set_alpha(0.0)
    fig.figurePatch.set_alpha(0.0)
    ax = fig.add_subplot(111)
    ax.patch.set_alpha(0.0)
    m=Basemap(projection='merc',llcrnrlat=26.5,urcrnrlat=lat[481,401-26],\
              llcrnrlon=-81.75,urcrnrlon=lon[481,401-26]+1,lat_ts=20,resolution='i')

    # Plot Data
    xi, yi = m(lon2,lat2)
    cs = m.pcolormesh(xi,yi,umag2,cmap=cm.gist_ncar_r,vmin=0.005, vmax = 0.15)
    xi, yi = m(lon,lat)
    cM = m.contour(   xi,yi,scope2,[1,1],colors='k',linewidth=2.0)
    
    # transform vectors to projection grid.
    uproj,vproj,xx,yy = m.transform_vector(u2[1::2,1::2]*7.5,v2[1::2,1::2]*7.5,\
                                           lonu[0,1::2].squeeze(),latu[1::2,0].squeeze(),\
                                           60,60,returnxy=True,masked=True)
                                           # now plot.
    Q  = m.quiver(xx,yy,uproj,vproj,scale=15)
    cm.gist_ncar.set_bad('white',alpha=0.0)

  # Save figure to file
  plot_file=plot_dir+'month_'+str(i+1)+'.png'
  #plot_file=plot_dir+'annual_mean.png'
  print 'Saving figure to file %s ' % plot_file
  plt.savefig(plot_file,dpi=150,bbox_inches='tight')

input_data.close()