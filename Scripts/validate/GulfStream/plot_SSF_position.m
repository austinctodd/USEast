%========================= plot_SSF_position.m ===========================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This reads the temperature field from the US East water age model    =%
%=  to find the location of the 10 degree isotherm at the surface (which =%
%=  is indicative of the Shelf/Slope Front).  The position is extracted  =%
%=  in the form of a latitude at each grid longitude (using interpolation=%
%=  routines.                                                            =%
%=                                                                       =%
%=========================================================================%

%--- Add various libraries and paths ---%
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/home/actodd/MYMATLAB/');
addpath(genpath('/home/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=========================================================================%
%=                  Set file input/output directories                    =%
%=========================================================================%
plotdir='/he_data/he/actodd/PLOTS/USeast/validate/KE/';
datadir='/gpfs_common/he_share/actodd/USeast-age/output/clim/';
matdir ='/he_data/he/actodd/DATA/GS_position/';

%=========================================================================%
%=                      Load ROMS grid information                       =%
%=========================================================================%
disp('Reading in ROMS output');

ncid=netcdf.open([datadir,'year1/useast_his_0001.nc'],'nowrite');
 lon        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'    ));
 lat        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'    ));
 mask       =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'   ));
 h          =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'          ));
netcdf.close(ncid);

%=========================================================================%
%=                        Load ROMS GS Positions                         =%
%=========================================================================%
lon2=lon(230:402,320);

load([matdir,'slope_front_lats.mat']);

tmeanlat=nanmean(templat,1);
tstdlat =nanstd(templat,1);

%=========================================================================%
%=                        Load ROMS GS Positions                         =%
%=========================================================================%
sdata=csvread([matdir,'SlopeFrontMonthlyMeansto1992.csv'],1,0);
sdata(find(sdata==0))=nan;

drink.lons=[-75:-50];
drink.slope_std =nanstd( sdata(:,3:2+length(drink.lons)));
drink.slope_mean=nanmean(sdata(:,3:2+length(drink.lons)));

%=========================================================================%
%=                    Plot grid and station locations                    =%
%=========================================================================%
figure(1); clf;
m_proj('lambert','long',[lon(230,320) lon(402,465)],...
                 'lat' ,[lat(230,320) lat(402,465)]);
m_pcolor(lon,lat,mask); shading flat; hold on;
caxis([-3 1]); colormap(gray);

%--- Bathymetry Contours ---%
m_contour(lon,lat,h.*mask,[100 500 1000 2000 3000 4000 ],'color',[.6 .6 .6],...
         'linewidth',1)
m_gshhs_i('line','color','k')                                     

%--- Site Locations ---%
h_Drink =m_patch([drink.lons drink.lons(end:-1:1) drink.lons(1)],...
                 [drink.slope_mean+drink.slope_std ...
		  drink.slope_mean(end:-1:1)-drink.slope_std(end:-1:1) ...
		  drink.slope_mean(1)+drink.slope_std(1)],[0.7 1 0.7],...
		 'edgecolor',[0 0.5 0]);%,'FaceAlpha',1);


hDrinkmean=m_plot(drink.lons,drink.slope_mean,'b','linewidth',2);
%          'markerfacecolor','b','markeredgecolor','k','markersize',5)

h_ROMS=m_patch([lon2' lon2(end:-1:1)' lon2(1)'],...
               [tmeanlat+tstdlat ...
		tmeanlat(end:-1:1)-tstdlat(end:-1:1) tmeanlat(1)+tstdlat(1)],...
		[1 0.7 0.7],'edgecolor',[.5 0 0],'FaceAlpha',0.95);
hROMSmean=m_plot(lon2',tmeanlat,'r','linewidth',2);
%hSTD_lat=m_plot(lon2',tmeanlat+tstdlat,'r--','linewidth',2);
%m_plot(lon2',tmeanlat-tstdlat,'r--','linewidth',2)

m_grid('tickdir','out');

ht=title('Mean Gulf Stream Position','fontsize',14,...
        'fontname','Helvetica','fontweight','bold');
aa=get(ht,'position');
set(ht,'position',[aa(1) aa(2)-0.007 aa(3)]);
set(gcf,'color','w');

%m_legend([hDrink,hMeanlat,hSTD_lat],'Drinkwater 1994','ROMS Mean position','ROMS Mean \pm \sigma','location','southeast');

%export_fig /he_data/he/actodd/PLOTS/USeast/GS_positions.png -png -r150
