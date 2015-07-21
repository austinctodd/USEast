%=============================================================================%
%=                                                                           =%
%=  Created by Austin C. Todd, NCSU, 2013                                    =%
%=                                                                           =%
%=  Program is designed to read in age data from the US East model output    =%
%=  and calculate the mean age over several different seasons. The mean age  =%
%=  will then be averaged over several time periods to correspond with calc- =%
%=  ulations done in Zhang et al., 2010 near the Hudson Bay shelf water      =%
%=  region. Averages will be for Spring (March - May), Summer (June-August), =%
%=  Fall (September - November), and Winter (December - February).           =%
%=                                                                           =%
%=============================================================================%

%--- Startup and setting of the path ---%
addpath('~/MYMATLAB/');
addpath(genpath('~/MYMATLAB/ROMS-matlab/'));
addpath(genpath('/Users/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=============================================================================%
%=                                                                           =%
%=       !!!!!!!!!!!        USER-DEFINED VARIABLES      !!!!!!!!!!           =%
%=                                                                           =%
%=============================================================================%
pdir='/Volumes/Black_box/Data/USeast-age/output/';

%=============================================================================%
%=                 Load in some grid and preliminary arrays                  =%
%=============================================================================%
refdate=datenum(2004,01,01); 

gridfile=['/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc'];

ncid=netcdf.open(gridfile,'nowrite');
h    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'         ));
lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'   ));
lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'   ));
mask1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_shelf'));
hmask=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'  ));
%time =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'));
netcdf.close(ncid)

lon=lon'; lat=lat';
nanmask1=mask1'./mask1';
nanmaskh=hmask'./hmask';

%========================================================================%
%=                     MABGOM Region - Multiple levels                      =%
%========================================================================%
ssnstrings=['JFM';'AMJ';'JAS';'OND'];
for ssn=1:4

  figure(1); clf;
%  ax1=axes('position',[0.1 0.05 0.8 0.93]);
  ax1=axes('position',[0.055 0.05 0.8 0.93]);
  m_proj('lambert','long',[-77 max(max(lon(1,1:end-25)))],...
                   'lat' ,[35 45.75]);
%m_contourf(lon,lat,hmask'); shading flat; colormap('gray'); caxis([-3 1])
  m_grid('tickdir','on'); hold on;
  m_pcolor(lon,lat,squeeze(seasonal_age3(:,:,end,ssn))'); 
  shading flat; caxis([50 365.25*2.5])
  [uu,vv]=regridromsvels2d(squeeze(seasonal_u(:,:,end,ssn))',...
                           squeeze(seasonal_v(:,:,end,ssn))');  
  m_quiver(lon(1:2:end,1:2:end),      lat(1:2:end,1:2:end),...
            uu(1:2:end,1:2:end).*nanmaskh(1:2:end,1:2:end),...
            vv(1:2:end,1:2:end).*nanmaskh(1:2:end,1:2:end),30,'k');
  %m_contour(lon,lat,hmask',[0 0],'k');
  m_contour(lon,lat,h',[50 100 250 500 1000 2000 3000 4000 5000],'k');
  %m_text(-81,44,'Mean suface vels','fontsize',14,'fontweight','bold');
  %m_text(-79.5,43.2,'(surface)','fontsize',14,'fontweight','bold');
  
  set(gcf,'color','w');
%wysiwyg;
 hcbar=colorbar('position',[0.85 0.05 0.1 0.9])
 set(hcbar,'position',[0.875 0.115 0.03 0.8])
 caxis([50 365.25*2.5]);
 set(hcbar,'linewidth',1,'fontsize',12,'fontweight','bold',...
     'ytick',[0.25:0.25:2.5].*365.25,'yticklabel',[0.25:0.25:2.5])
 set(get(hcbar,'ylabel'),'String','Mean Age (yr)','fontsize',14)


  print('-dpng','-r150','-painters',['/Volumes/Black_box/Data/PLOTS/',...
         'USeast-age/vels/',ssnstrings(ssn,:),'_age_MABGOM.png']);
     pause(0.5)
end

%========================================================================%
%=                  Gulf of Mexico - Multiple levels                    =%
%========================================================================%
ssnstrings=['JFM';'AMJ';'JAS';'OND'];
for ssn=1:4

  figure(1); clf;
  ax1=axes('position',[0.055 0.05 0.8 0.93]);
  m_proj('lambert','long',[-98 -79],'lat' ,[17 32]);
%m_contourf(lon,lat,hmask'); shading flat; colormap('gray'); caxis([-3 1])
  m_grid('tickdir','on'); hold on;
  m_pcolor(lon,lat,squeeze(seasonal_age2(:,:,end,ssn))'); 
  shading flat; caxis([30 365.25*2])
  [uu,vv]=regridromsvels2d(squeeze(seasonal_u(:,:,end,ssn))',...
                           squeeze(seasonal_v(:,:,end,ssn))');  
  m_quiver(lon(1:2:end,1:2:end),      lat(1:2:end,1:2:end),...
            uu(1:2:end,1:2:end).*nanmaskh(1:2:end,1:2:end),...
            vv(1:2:end,1:2:end).*nanmaskh(1:2:end,1:2:end),30,'k');
  %m_contour(lon,lat,hmask',[0 0],'k');
  m_contour(lon,lat,h',[50 100 250 500 1000 2000 3000 4000 5000],'k');
  %m_text(-81,44,'Mean suface vels','fontsize',14,'fontweight','bold');
  %m_text(-79.5,43.2,'(surface)','fontsize',14,'fontweight','bold');
  
  set(gcf,'color','w');
 % wysiwyg;
 
 hcbar=colorbar('position',[0.85 0.05 0.1 0.9])
 set(hcbar,'position',[0.875 0.115 0.03 0.8])
 caxis([30 365.25*2]);
 set(hcbar,'linewidth',1,'fontsize',12,'fontweight','bold',...
     'ytick',[0.25:0.25:2].*365.25,'yticklabel',[0.25:0.25:2])
 set(get(hcbar,'ylabel'),'String','Mean Age (yr)','fontsize',14)
 
 
  print('-dpng','-r150','-painters',['/Volumes/Black_box/Data/PLOTS/',...
         'USeast-age/vels/',ssnstrings(ssn,:),'_age2_GOM.png']);
     pause(0.5)
end




