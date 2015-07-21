#========================= plot_15C_isotherm.m ===========================#
#=                                                                       =#
#=  Written by Austin C Todd, NCSU (2014)                                =#
#=                                                                       =#
#=  This reads the temperature field from the US East water age model    =#
#=  to find the location of the 15 degree isotherm at 200m depth (which  =#
#=  is indicative of the Gulf Stream Wall).  The position is extracted   =#
#=  in the form of a latitude at each grid longitude (using interpolation=#
#=  routines.                                                            =#
#=                                                                       =#
#=========================================================================#

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


#=========================================================================#
#=                  Set file input/output directories                    =#
#=========================================================================#
plotdir    = '/Volumes/Black_box/Data/PLOTS/USeast-age/validate/';
ROMSGSfile = '/Volumes/Black_box/Data/USeast/Data/GS_position/'+\
             'ROMS_GS_lats.txt';
Drinkwater = '/Volumes/Black_box/Data/USeast/Data/GS_position/'+\
             'FrontalMonthlyMeansto1992.csv';
grid_file  = '/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc'

#=========================================================================#
#=                      Load ROMS grid information                       =#
#=========================================================================#
print 'Ingesting data from file %s ...' % (grid_file)

grid_data=Dataset(grid_file,mode='r')
mask= grid_data.variables['mask_rho' ][:,0:401-25]
h    = grid_data.variables['h'        ][:,0:401-25]
lon =grid_data.variables['lon_rho'    ][:,0:401-25]
lat =grid_data.variables['lat_rho'    ][:,0:401-25]
grid_data.close()

#=========================================================================#
#=                        Load ROMS GS Positions                         =#
#=========================================================================#
lon2=lon[289,201:402]

f=open(ROMSGSfile,'r')

for line in f:
  line=line.strop()
  columns=line.split()
  print(float(columns[0]),float(columns[1])
#=========================================================================#
#=                        Load ROMS GS Positions                         =#
#=========================================================================#
drinkwater=textread([matdir,'drinkwater1994.txt']);
data=csvread([matdir,'FrontalMonthlyMeansto1992.csv'],1,0);
data(find(data==0))=nan;

drink.lons=[-75:-50];
drink.std =nanstd( data(:,3:2+length(drink.lons)));
drink.mean=nanmean(data(:,3:2+length(drink.lons)));

sdata=csvread([matdir,'SlopeFrontMonthlyMeansto1992.csv'],1,0);
sdata(find(sdata==0))=nan;

drink.slope_std =nanstd( sdata(:,3:2+length(drink.lons)));
drink.slope_mean=nanmean(sdata(:,3:2+length(drink.lons)));

#=========================================================================#
#=                    Plot grid and station locations                    =#
#=========================================================================#
figure(1); clf;
m_proj('lambert','long',[-78 max(lon(:))],...
                 'lat' ,[33 43]);
m_pcolor(lon,lat,mask); shading flat; hold on;
caxis([-3 1]); colormap(gray);

#--- Bathymetry Contours ---#
m_contour(lon,lat,h.*mask,[100 500 1000 2000 3000 4000 ],'color',[.6 .6 .6],...
         'linewidth',1)
m_gshhs_i('line','color','k')                                     

#--- Site Locations ---#
h_Drink=m_patch([drink.lons drink.lons(end:-1:1) drink.lons(1)],...
                [drink.mean+drink.std ...
		 drink.mean(end:-1:1)-drink.std(end:-1:1) drink.mean(1)+drink.std(1)],[0.7 0.7 1],...
		 'edgecolor',[0 0 .5]);#,'FaceAlpha',1);


hDrinkmean=m_plot(-drinkwater(:,1),mean(drinkwater(:,2:end),2),'b','linewidth',2);
#          'markerfacecolor','b','markeredgecolor','k','markersize',5)

h_ROMS=m_patch([lon2(25:end)' lon2(end:-1:25)' lon2(25)'],...
                [tmeanlat(25:end)+tstdlat(25:end) ...
		 tmeanlat(end:-1:25)-tstdlat(end:-1:25) tmeanlat(25)+tstdlat(25)],[1 0.7 0.7],...
		 'edgecolor',[.5 0 0],'FaceAlpha',0.95);
hROMSmean=m_plot(lon2',tmeanlat,'r','linewidth',2);
#hSTD_lat=m_plot(lon2',tmeanlat+tstdlat,'r--','linewidth',2);
#m_plot(lon2',tmeanlat-tstdlat,'r--','linewidth',2)

m_grid('tickdir','out');

ht=title('Mean Gulf Stream Position','fontsize',14,...
        'fontname','Helvetica','fontweight','bold');
aa=get(ht,'position');
set(ht,'position',[aa(1) aa(2)-0.007 aa(3)]);
set(gcf,'color','w');

#m_legend([hDrink,hMeanlat,hSTD_lat],'Drinkwater 1994','ROMS Mean position','ROMS Mean \pm \sigma','location','southeast');

export_fig /gpfs_backup/he_data/he/actodd/PLOTS/USeast/GS_position_new.png -png -r150
