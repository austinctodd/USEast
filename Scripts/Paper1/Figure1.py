#===============================================================================
#
# SCRIPT :  plot_domains.py
#
# PURPOSE : Ingest grid data from US East water age model and plots the model
#           domain on a map with various other lines and labels.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 1 JULY 2015.
#
#===============================================================================

input_file = '/Volumes/Black_box/Data/USeast-age/output/clim/averages/avg_3hrly.nc'
grid_file  = '/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc'
ATL_file   = '/Volumes/Black_box/Data/USeast/Data/grd/grid_ATLscope.nc'
GOM_file   = '/Volumes/Black_box/Data/USeast/Data/grd/grid_GOM_shelf_scope.nc'
cmap_file  = '/Users/actodd/Documents/ROMS-utils/USeast/ElevationBathymetryColormap/etopo1.clr'
river_file = '/Volumes/Black_box/Data/USeast/Data/frc/USeast-riverNOVN-12.nc'
plot_file  = '/Volumes/Black_box/Data/PLOTS/USeast-age/domain/domain.png'

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
print 'Ingesting data from file %s ...' % (grid_file)
grid_data=Dataset(grid_file,mode='r')
hmask= grid_data.variables['mask_rho' ][:,0:401-25]
Mmask= grid_data.variables['mask_Miss'][:,0:401-25]
h    = grid_data.variables['h'        ][:,0:401-25]
hraw = grid_data.variables['hraw'     ][0,:,0:401-25]
lon =grid_data.variables['lon_rho'    ][:,0:401-25]
lonu=grid_data.variables['lon_u'      ][:]
lonv=grid_data.variables['lon_v'      ][:]
lat =grid_data.variables['lat_rho'    ][:,0:401-25]
latu=grid_data.variables['lat_u'      ][:]
latv=grid_data.variables['lat_v'      ][:]
grid_data.close()

#-------------------------------------------------------------------------------
# Open ATL/GOM/GOM shelf files and read in scopes
#-------------------------------------------------------------------------------
print 'Ingesting data from file %s ...' % (ATL_file)
grid_data=Dataset(ATL_file,mode='r')
Ascope= grid_data.variables['scope_rho' ][:,0:401-25]
grid_data.close()

print 'Ingesting data from file %s ...' % (GOM_file)
grid_data=Dataset(GOM_file,mode='r')
Gscope= grid_data.variables['scope_rho' ][:,0:401-25]
grid_data.close()

Sscope=Ascope+Gscope

#-------------------------------------------------------------------------------
# Open river file and read in variables
#-------------------------------------------------------------------------------
print 'Ingesting data from file %s ...' % (river_file)
river_data=Dataset(river_file,mode='r')
river_X    =river_data.variables['river_Xposition'][:]
river_Y    =river_data.variables['river_Eposition'][:]
river_trans=river_data.variables['river_transport'][:]
river_dir  =river_data.variables['river_direction'][:]
river_flag =river_data.variables['river_flag'     ][:]
river_data.close()

river_lon =np.zeros(shape=(len(river_X)))
river_lat =np.zeros(shape=(len(river_X)))
river_size=np.zeros(shape=(len(river_X)))
river_trans=np.mean(river_trans,axis=0)

#-------------------------------------------------------------------------------
# Set river mouth lat/lons and scale direction
#-------------------------------------------------------------------------------
for i in range(0,len(river_X)):
  if river_dir[i]:
    river_lon[ i]=lonv[river_Y[i],river_X[i]]
    river_lat[ i]=latv[river_Y[i],river_X[i]]
    river_size[i]=np.abs(river_trans[i])
  else:
    river_lon[ i]=lonu[river_Y[i],river_X[i]]
    river_lat[ i]=latu[river_Y[i],river_X[i]]
    river_size[i]=np.abs(river_trans[i])


#-------------------------------------------------------------------------------
# Read in colormap for bathymetry
#-------------------------------------------------------------------------------
#bathymap=np.loadtxt(cmap_file,skiprows=2,usecols=(1,2,3))
#bathymap=bathymap/255
#cmap1 = cm.colors.LinearSegmentedColormap.from_list("my_colormap",bathymap[0:11999,:])

#-------------------------------------------------------------------------------
# Mask the age variables to the shelf
#-------------------------------------------------------------------------------
hraw=np.where(hraw>0.0,0.0,1.0)
for i in range(0,482):
  for j in range(0,401-25):
    if (hraw[i,j]==1.0):
      if (hmask[i,j]==1.0):
        hraw[i,j]=0.0
#hraw=1-hraw
h=np.ma.array(h,mask=hraw)
#age=age/86400.
cm.ocean.set_bad('white',alpha=0)

Gmask=1-Gscope
Gmask=np.ma.array(Gmask,mask=Gmask)

Smask=1-Sscope
Smask=np.ma.array(Smask,mask=Smask)

#-------------------------------------------------------------------------------
# Plot Data using Basemap mercator projection
#-------------------------------------------------------------------------------
print 'Plotting data.'

fig=plt.figure(num=None, figsize=(10, 10), dpi=150, facecolor='w')
m=Basemap(projection='merc',llcrnrlat=lat.min(),urcrnrlat=lat.max(),\
          llcrnrlon=lon.min(),urcrnrlon=lon.max(),lat_ts=20,resolution='i')

# Draw a thick border around the whole map
m.drawmapboundary(fill_color='white',linewidth=3)

# Plot Data
xi, yi = m(lon,lat)
cs = m.pcolormesh(xi,yi,-h,cmap=cm.ocean,vmin=-7000, vmax=150)
Mc1= m.contourf(xi,yi,Smask,levels=[0,0],cmap=cm.gnuplot,vmin=-1, vmax=0, alpha=0.5)
#Mc = m.contourf(xi,yi,Smask,levels=[0,0],cmap=cm.jet,vmin=-1, vmax=0)
cM1 = m.contour(   xi,yi,Sscope,[0.5,0.5],colors=[(0.6,0.6,0)],linewidth=2)
#cM = m.contour(   xi,yi,Sscope,[0.5,0.5],colors='k',linewidth=2.0)

# Contour depths (but first mask Pacific)
h=np.ma.array(h,mask=1.0-hmask)
cn = m.contour(xi,yi,h,[100,500,1000,2000,3000,4000],linewidth=1.5,colors='0.7')

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

# Add river points scattered and sized by river output
x,y = m(river_lon,river_lat)
for j in range(0,len(river_lon)):
    m.plot(x[j],y[j],markersize=np.log(river_size[j])*3,marker='o',color='r',alpha=0.4)

# Add Rectangle as base for legend
lns=[-98.5,-87.25,-87.25,-98.5]
#lts=[38,38,46.1,46.1]
lts=[35,35,46.1,46.1]
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

# Print text next to each circle
lns=[-94.25,-94.25,-94.25,-94.25,-94.25,-94.25]
lts=[ 44.35, 43.35, 42.25, 40.9, 39.7, 38.85]
x1,y1 = m(lns,lts)
plt.text(x1[0],y1[0],'5 m$^3$/s',fontsize=12,ha='left',va='center',color='k')
plt.text(x1[1],y1[1],'50 m$^3$/s',fontsize=12,ha='left',va='center',color='k')
plt.text(x1[2],y1[2],'500 m$^3$/s',fontsize=12,ha='left',va='center',color='k')
plt.text(x1[3],y1[3],'5000 m$^3$/s',fontsize=12,ha='left',va='center',color='k')

# Print legend title
x1,y1 = m([-92.875],[45.25])
plt.text(x1[0],y1[0],'Mean Streamflow',fontsize=12,fontweight='bold',ha='center',va='center',color='k')

# Create patch as a bathymetry colorbar
a  =np.arange(0,1.01,0.01)
lns=np.array([[(a*(-88.75+97.25))-97.25],[(a*(-88.75+97.25))-97.25]]).squeeze()
lts=np.array([[a*0+36.25],[a*0+38.25]]).squeeze()
bth=np.array([[a*-4000],[a*-4000]]).squeeze()
x2,y2 = m(lns,lts)
cbs=m.pcolormesh(x2,y2,bth, cmap=cm.ocean, vmin=-7000, vmax=150)
cbs.set_zorder(4)

# Make lines around Colorbar
m.plot(x2[0,:],y2[0,:],color='k',linewidth=1,zorder=4)
m.plot(x2[1,:],y2[1,:],color='k',linewidth=1,zorder=4)
m.plot(x2[:,0],y2[:,0],color='k',linewidth=1,zorder=4)
m.plot(x2[:,100],y2[:,100],color='k',linewidth=1,zorder=4)

# Write ticks below colorbar
plt.text(x2[0,0],y2[0,0]-15000,'0',fontsize=10,ha='center',va='top',color='k')
plt.text(x2[0,25],y2[0,25]-15000,'1',fontsize=10,ha='center',va='top',color='k')
plt.text(x2[0,50],y2[0,50]-15000,'2',fontsize=10,ha='center',va='top',color='k')
plt.text(x2[0,75],y2[0,75]-15000,'3',fontsize=10,ha='center',va='top',color='k')
plt.text(x2[0,100],y2[0,100]-15000,'4',fontsize=10,ha='center',va='top',color='k')


# Print Colorbar title
x1,y1 = m([-92.875,-92.875],[38.9,38.05])
plt.text(x1[0],y1[0],'Model Depth (km)',fontsize=12,fontweight='bold',ha='center',va='center',color='k')
#plt.text(x1[1],y1[1],'(km)',fontsize=12,fontweight='bold',ha='center',va='center',color='k')

# Add Colorbar
#clabels=np.array([0,0.25,0.5,0.75,1,1.25,1.5,1.75,2])
#clabels=clabels*365.25
#cbar = m.colorbar(cs, location='bottom', pad="4%")#,ticks=clabels)
#cbar.set_label('Mean Age (years)')
#cbar.ax.set_xticklabels([0,0.25,0.5,0.75,1,1.25,1.5,1.75,2])

# Add Title
#plt.title('Age tracer 3')

# Save figure to file
print 'Saving figure to file %s ' % plot_file
plt.savefig(plot_file,dpi=150,bbox_inches='tight')

