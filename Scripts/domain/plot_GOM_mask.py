#===============================================================================
#
# SCRIPT :  plot_domains.py
#
# PURPOSE : Ingest grid data from US East water age model and plots the model
#           domain on a map with various other lines and labels.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 3 FEBRUARY 2015.
#
#===============================================================================

input_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/avg_3hrly.nc'
GOM_grid_file  = '/Volumes/Black_box/Data/USeast/Data/grd/grid_ATLscope.nc'
GOM_shelf_file  = '/Volumes/Black_box/Data/USeast/Data/grd/grid_ATLscope.nc'
cmap_file  = '/Users/actodd/Documents/ROMS-utils/USeast/ElevationBathymetryColormap/etopo1.clr'
river_file  = '/Volumes/Black_box/Data/USeast/Data/frc/USeast-riverNOVN-12.nc'
plot_file  = '/Volumes/Black_box/Data/PLOTS/USeast-age/domain/map_ATL.png'

#-------------------------------------------------------------------------------
# Define all required libraries, routines, and modules
#-------------------------------------------------------------------------------
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import sys
import mpl_toolkits.basemap
import os
import PIL

#-------------------------------------------------------------------------------
# Define all required sub-libraries, sub-routines, and sub-modules
#-------------------------------------------------------------------------------
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
from matplotlib.patches import Polygon

#-------------------------------------------------------------------------------
# Print message to user
#-------------------------------------------------------------------------------
print 'Ingesting data from file %s ...' % (input_file)

#-------------------------------------------------------------------------------
# Open files and read in variables
#-------------------------------------------------------------------------------
GOM_grid_data=Dataset(GOM_grid_file,mode='r')
hmask= GOM_grid_data.variables['mask_rho'  ][:,0:401-25]
Mmask= GOM_grid_data.variables['scope_rho' ][:,0:401-25]
h    = GOM_grid_data.variables['h'         ][:,0:401-25]
hraw = GOM_grid_data.variables['hraw'      ][0,:,0:401-25]
lon  = GOM_grid_data.variables['lon_rho'   ][:,0:401-25]
lonu = GOM_grid_data.variables['lon_u'     ][:]
lonv = GOM_grid_data.variables['lon_v'     ][:]
lat  = GOM_grid_data.variables['lat_rho'   ][:,0:401-25]
latu = GOM_grid_data.variables['lat_u'     ][:]
latv = GOM_grid_data.variables['lat_v'     ][:]
GOM_grid_data.close()

GOM_shelf_data=Dataset(GOM_shelf_file,mode='r')
Smask= GOM_shelf_data.variables['scope_rho'][:,0:401-25]
GOM_shelf_data.close()

#-------------------------------------------------------------------------------
# Mask the age variables to the shelf
#-------------------------------------------------------------------------------
#Mmask=Mmask*(1-Smask)
Mmask=1-Mmask
hmask=1-hmask
Smask=1-Smask
h=np.ma.array(h,mask=hmask)
Mmask=np.ma.array(Mmask,mask=Mmask)
Smask=np.ma.array(Smask,mask=Smask)
#age=age/86400.
#cm.ocean.set_bad('white',alpha=0)
#cm.jet.set_bad('white',alpha=0)

#-------------------------------------------------------------------------------
# Plot Data using Basemap mercator projection
#-------------------------------------------------------------------------------
print 'Plotting data.'

fig=plt.figure(num=None, figsize=(10, 10), dpi=150, facecolor='w')
m=Basemap(projection='merc',llcrnrlat=17,urcrnrlat=31.5,\
          llcrnrlon=-98.5,urcrnrlon=-79,lat_ts=20,resolution='f')
m=Basemap(projection='merc',llcrnrlat=27,urcrnrlat=lat.max(),\
          llcrnrlon=-82,urcrnrlon=lon.max(),lat_ts=20,resolution='i')

# Draw a thick border around the whole map
m.drawmapboundary(fill_color='white',linewidth=3)

# Plot Data
xi, yi = m(lon,lat)
cs = m.pcolormesh(xi,yi,-h,cmap=cm.ocean,vmin=-7000, vmax=150)

# Contour depths (but first mask Pacific)
cn = m.contour(xi,yi,h,[100,500,1000,2000,3000,4000],linewidth=1.5,colors='0.7')
#Mm = m.contourf(xi,yi,1-Mmask,levels=[1,1],cmap=cm.cool,vmin=0, vmax=1,alpha=0.6)
Mc = m.contourf(xi,yi,Smask,levels=[0,0],cmap=cm.jet,vmin=-1, vmax=0, alpha=0.6)
Mc2= m.contour( xi,yi,Smask,levels=[0,0],colors='r',linewitdh=2)
#Mm2= m.contour( xi,yi,1-Mmask,levels=[0,0],colors='0.2',linewitdh=2)

# Add Grid Lines
m.drawparallels(np.arange(0.,90.,5.), labels=[1,0,0,0],\
                fontsize=10, linewidth=0.75, color='.5')
m.drawmeridians(np.arange(-105.,-35.,5.), labels=[0,0,0,1],\
                fontsize=10, linewidth=0.75, color='.5')

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.etopo()
#m.bluemarble()
#m.fillcontinents(color='grey')

# Add river points scattered and sized by river output
#x,y = m(river_lon,river_lat)
  #for j in range(0,len(river_lon)):
#  if river_flag[j]>7:
#    m.plot(x[j],y[j],markersize=np.log(river_size[j])*2,marker='o',color='y',alpha=0.4)


# Save figure to file
print 'Saving figure to file %s ' % plot_file
plt.savefig(plot_file,dpi=150,bbox_inches='tight')

