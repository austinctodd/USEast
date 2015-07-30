%============================ plot_site_locations.m ===========================%
%=                                                                            =%
%=  Written by Austin C. Todd, NCSU (2014) for personal use                   =%
%=                                                                            =%
%=  Program loads US East Coast ROMS Bathymetry and grid and plots the site   =%
%=  locations that will be used for particle advection experiments.           =%
%=                                                                            =%
%==============================================================================%

addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/home/actodd/MYMATLAB/');
addpath(genpath('/home/actodd/MYMATLAB/m_map'));
format long g; format compact;

%==============================================================================%
%=                     Set file input/output directories                      =%
%==============================================================================%
plotdir  ='/he_data/he/actodd/PLOTS/Connectivity/';
romsdir  ='/gpfs_common/he_share/actodd/USeast-age/output/clim/';
ltransdir='/home/actodd/LTRANS/output/site';

%==============================================================================%
%=                         Load ROMS grid information                         =%
%==============================================================================%
disp('Reading in ROMS output');

ncid=netcdf.open([romsdir,'year1/useast_his_0001.nc'],'nowrite');
 lon        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'    ));
 lat        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'    ));
 mask       =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'   ));
 h          =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'          ));
netcdf.close(ncid);

%==============================================================================%
%=       Set the site lat/lon pairs (NOTE: sites 4-6 will not be used)        =%
%==============================================================================%
sitelat=[26.354722, 27.723056, 26.03,      11.233367, 10.327833,...
	 10.327833, 32.490833, 32.9741617, 37.0985,   38.166667];
sitelon=[-94.496667, -91.275,    -84.915,    -59.34595,  -58.88895, ...
         -58.88895,  -76.188333, -75.921667, -74.623167, -73.833333];
sitenames={'Alaminos Canyon';'Brine Pool';'Florida Escarpment';'El Pilar';...
           'Orenoque A';'Orenoque B ';'Blake Ridge';'Cape Fear';...
	   'Norfolk Canyon';'Baltimore Canyon'};

sitelat=sitelat([1:3 7:10]); sitelon=sitelon([1:3 7:10]); 
sitenames=sitenames([1:3 7:10]);

sitelat(1)=lat( 47,210); sitelon(1)=lon( 47,210);
sitelat(3)=lat(144,201); sitelon(3)=lon(144,201);

%--- Set colors of each site ---%
cols=[0.7 0 0; 0.7 0.4 0; 0.4 0.2 0.6; 0 0 0.8; 0 0.4 0; 0.6 0.3 0; 1 0 0.5];

%==============================================================================%
%=                           Read in LTRANS output                            =%
%==============================================================================%
%ncid=netcdf.open([ltransdir,'site3.nc'],'nowrite');
%  ltrans.lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
%  ltrans.lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
%  ltrans.dob  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dob'));
%  ltrans.age  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age'));
%  ltrans.depth=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'depth'));
%  ltrans.color=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'color'));
%netcdf.close(ncid);

%==============================================================================%
%=                           Arrange output by dob                            =%
%==============================================================================%


%==============================================================================%
%=                       Plot grid and station locations                      =%
%==============================================================================%
for site=2:2
  ncid=netcdf.open([ltransdir,sprintf('%1i',site),'/site',sprintf('%1i',site),'.nc'],'nowrite');
    ltrans.lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
    ltrans.lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
  netcdf.close(ncid);

  figure(1); clf;
  ax1=axes('position',[0.1 0.05 0.8 0.93]);
%  m_proj('lambert','long',[min(lon(:)) max(max(lon(1:end-25,:)))],...
%                   'lat' ,[min(lat(:)) max(max(lat(1:end-25,:)))]);
  m_proj('lambert','long',[-98 -89],...
                    'lat' ,[25 31]);
  m_pcolor(lon,lat,mask); shading flat; hold on;
  colormap(gray);
  caxis([-4 1]);

  %--- Bathymetry Contours ---%
  m_contour(lon,lat,h.*mask,[100 500 1000 2000 3000],'color',[.6 .6 .6],...
           'linewidth',1)
  m_contour(lon,lat,mask,[0 0],'k','linewidth',1)                                     


  for i=1:365
     m_plot(ltrans.lon(i,i:i+27),ltrans.lat(i,i:i+27),'color',cols(site,:))
%   m_plot(ltrans.lon(i,i+28:i+181),ltrans.lat(i,i+28:i+181),'color',[0 .4 0])
%   m_plot([sitelon(2) ltrans.lon(i,i:i+27)],[sitelat(2) ltrans.lat(i,i:i+27)],'b')
  end

  %--- Site Locations ---%
  m_plot(sitelon(site),sitelat(site),'rs','markersize',5,'markerfacecolor',...
         'r','markeredgecolor','k')
  m_text(-97.5,30,sitenames{site},'color',cols(site,:),'fontsize',12,...
         'fontname','helvetica','fontweight','bold');
  m_text(-97.5,29.6,'4 Weeks advection','color',cols(site,:),'fontsize',12,...
         'fontname','helvetica','fontweight','bold');

  m_grid('tickdir','out')
  set(gcf,'color','w');
  eval(['export_fig /he_data/he/actodd/PLOTS/Connectivity/site_',...
        sprintf('%1i',site),'_4weeks.png -png -r150 -painters']);
end

return;
%--- Loop through to put site numbers on map ---%
for i=1:length(sitelon)
  if (i==5||i==3)
  m_text(sitelon(i)+0.75,sitelat(i),sprintf('%1i',i),'color','r',...
          'fontsize',12,'fontweight','bold','fontname','helvetica');
  else
  m_text(sitelon(i)+0.75,sitelat(i)-0.75,sprintf('%1i',i),'color','r',...
          'fontsize',12,'fontweight','bold','fontname','helvetica');
  end
end

m_grid('tickdir','out')

ax2=axes('position',[0.68 0.075 0.3 0.9]);
set(ax2,'position',[0.145 0.23 0.3 0.9]);
%plot(ones(1,7).*-1.25,[7.25:-.5:4.25],'rs','markersize',10,'markeredgecolor','k',...
%                       'markerfacecolor','r');
patch([-1.45 1.58 1.58 -1.45],[4.4 4.4 8 8],'w');
hold on;
for i=1:7	       
  text(-1.37,7.73-(i*0.44),[sprintf('%1i',i),'. ',sitenames{i}],'color','r',...
        'fontsize',12,'fontname','helvetica','fontweight','bold');
end
axis([-2 2 0 10]); axis off
text(-0.2,7.73,'Site Locations','fontsize',14,'fontweight','bold',...
            'fontname','Helvetica','HorizontalAlignment','Center');
%export_fig /he_data/he/actodd/PLOTS/Connectivity/site_locations.png -png -r150 -painters
