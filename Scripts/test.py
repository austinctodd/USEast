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

grid_file  = '/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc'
plot_file  = '/Volumes/Black_box/Data/PLOTS/USeast-age/domain/temp.png'

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
# Open files and read in variables
#-------------------------------------------------------------------------------
grid_data=Dataset(grid_file,mode='r')
lon =grid_data.variables['lon_rho'    ][:,0:401-25]
lat =grid_data.variables['lat_rho'    ][:,0:401-25]
grid_data.close()

#-------------------------------------------------------------------------------
# Plot Data using Basemap mercator projection
#-------------------------------------------------------------------------------
print 'Plotting data.'

fig=plt.figure(num=None, figsize=(10, 10), dpi=150, facecolor='w')
m=Basemap(projection='merc',llcrnrlat=lat.min(),urcrnrlat=lat.max(),\
          llcrnrlon=lon.min(),urcrnrlon=lon.max(),lat_ts=20,resolution='i')

# Draw a thick border around the whole map
m.drawmapboundary(fill_color='white',linewidth=3)

# Add Grid Lines
m.drawparallels(np.arange(0.,90.,10.), labels=[1,0,0,0],\
                fontsize=10, linewidth=0.75, color='.5')
m.drawmeridians(np.arange(-105.,-35.,10.), labels=[0,0,0,1],\
                fontsize=10, linewidth=0.75, color='.5')

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.etopo()
#m.bluemarble()
#m.fillcontinents(color='grey')

# Add Rectangle as base for legend
lns=[-98.5,-87.25,-87.25,-98.5]
#lts=[38,38,46.1,46.1]
lts=[33,33,46.1,46.1]
x,y = m(lns,lts)
xy = zip(x,y)
poly = Polygon( xy, facecolor='white')
pt= plt.gca().add_patch(poly)
pt.set_zorder(3)

# Add Scale circles for legend
lns=[-96,-96,-96,-96,-96]
lts=[ 44.35, 43.45, 42.35, 40.9, 39.25]
x1,y1 = m(lns,lts)
m.plot(x1[0],y1[0],markersize=np.log(   5.0)*3,marker='o',color='r',alpha=0.6,zorder=4)
m.plot(x1[1],y1[1],markersize=np.log(  50.0)*3,marker='o',color='r',alpha=0.6,zorder=4)
m.plot(x1[2],y1[2],markersize=np.log( 500.0)*3,marker='o',color='r',alpha=0.6,zorder=4)
m.plot(x1[3],y1[3],markersize=np.log(5000.0)*3,marker='o',color='r',alpha=0.6,zorder=4)
m.plot(x1[4],y1[4],markersize=np.log(5000.0)*3,marker='o',color='y',alpha=0.6,zorder=4)

# Print text next to each circle
lns=[-94.25,-94.25,-94.25,-94.25,-94.25,-94.25]
lts=[ 44.35, 43.35, 42.25, 40.9, 39.7, 38.85]
x1,y1 = m(lns,lts)
plt.text(x1[0],y1[0],'5 m$^3$/s',fontsize=12,ha='left',va='center',color='k')
plt.text(x1[1],y1[1],'50 m$^3$/s',fontsize=12,ha='left',va='center',color='k')
plt.text(x1[2],y1[2],'500 m$^3$/s',fontsize=12,ha='left',va='center',color='k')
plt.text(x1[3],y1[3],'5000 m$^3$/s',fontsize=12,ha='left',va='center',color='k')
plt.text(x1[4],y1[4],'Mississippi /',fontsize=12,ha='left',va='center',color='k')
plt.text(x1[5],y1[5],'Atchafalaya',fontsize=12,ha='left',va='center',color='k')

# Print legend title
x1,y1 = m([-92.875],[45.25])
plt.text(x1[0],y1[0],'Mean Streamflow',fontsize=12,fontweight='bold',ha='center',va='center',color='k')

# Create patch as a bathymetry colorbar
a  =np.arange(0,1.01,0.01)
lns=np.array([[(a*(-88.75+97.25))-97.25],[(a*(-88.75+97.25))-97.25]]).squeeze()
lts=np.array([[a*0+34],[a*0+36]]).squeeze()
bth=np.array([[(a*4000)-4000],[(a*4000)-4000]]).squeeze()
x2,y2 = m(lns,lts)
cbs=m.pcolormesh(x2,y2,bth, cmap=cm.ocean, vmin=-7000, vmax=150)
cbs.set_zorder(4)

# Make lines around Colorbar
m.plot(x2[0,:],y2[0,:],color='k',linewidth=1,zorder=4)
m.plot(x2[1,:],y2[1,:],color='k',linewidth=1,zorder=4)
m.plot(x2[:,0],y2[:,0],color='k',linewidth=1,zorder=4)
m.plot(x2[:,100],y2[:,100],color='k',linewidth=1,zorder=4)

# Write ticks below colorbar
plt.text(x2[0,0],y2[0,0]-8000,'4',fontsize=10,ha='center',va='top',color='k')
plt.text(x2[0,25],y2[0,25]-8000,'3',fontsize=10,ha='center',va='top',color='k')
plt.text(x2[0,50],y2[0,50]-8000,'2',fontsize=10,ha='center',va='top',color='k')
plt.text(x2[0,75],y2[0,75]-8000,'1',fontsize=10,ha='center',va='top',color='k')
plt.text(x2[0,100],y2[0,100]-8000,'0',fontsize=10,ha='center',va='top',color='k')


# Print Colorbar title
x1,y1 = m([-92.875,-92.875],[37.4,36.55])
plt.text(x1[0],y1[0],'Model Bathymetry',fontsize=12,fontweight='bold',ha='center',va='center',color='k')
plt.text(x1[1],y1[1],'(km)',fontsize=12,fontweight='bold',ha='center',va='center',color='k')

# Save figure to file
print 'Saving figure to file %s ' % plot_file
plt.savefig(plot_file,bbox_inches='tight')

