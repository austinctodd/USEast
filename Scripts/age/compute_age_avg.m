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
%=  For each source region, there are two tracer types:                      =%
%=             C = water age tracer concentration                            =%
%=         alpha = age concentration                                         =%
%=                                                                           =%
%=  where                                                                    =%
%=                                                                           =%
%=         dC/dt = P - D - del*(UC       - K*(delC      ))                   =%
%=   d(alpha)/dt =     C - del*(U(alpha) - K*(del(alpha)))                   =%
%=                                                                           =%
%=  P = source, D = sink, del = d/dx+d/dy, * = dot product, U = velocity,    =%
%=  and K = diffusivity                                                      =%
%=                                                                           =%
%=                      The mean age A = alpha / C                           =%
%=                                                                           =%
%=============================================================================%

%--- Startup and setting of the path ---%
addpath('~/MYMATLAB/');
addpath(genpath('~/MYMATLAB/ROMS-matlab/'));
addpath(genpath('/home/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=============================================================================%
%=                                                                           =%
%=       !!!!!!!!!!!        USER-DEFINED VARIABLES      !!!!!!!!!!           =%
%=                                                                           =%
%=============================================================================%
tracer='shelf'; % 'shelf' or 'deep' or 'river'
pdir='/gpfs_share/actodd/USeast-age/output/2005';

%=============================================================================%
%=                  Setup some grid and preliminary arrays                   =%
%=============================================================================%
refdate=datenum(2004,01,01); % Why the hell do we choose this date???

gridfile=[pdir,'/useast_his.nc'];

ncid=netcdf.open(gridfile,'nowrite');
h    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'         ));
lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'   ));
lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'   ));
mask1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_shelf'));
hmask=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'  ));
time =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'));
%netcdf.close(ncid)

lon=lon'; lat=lat';
nanmask1=mask1'./mask1';
nanmaskh=hmask'./hmask';
time=time./86400+refdate;

%=============================================================================%
%=           Loop to read in files and calculate surface mean age            =%
%=============================================================================%
refdate=datenum(2004,01,01); 

%--- Spring 2005 ---%
stind=find(time==datenum(2005,3,1))-1;
enind=find(time==datenum(2005,6,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_06'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
mean_age=squeeze(double(var2./var1));
mean_age(find(var1<1e-6))=nan;

spring.ma(1,:,:)=squeeze(nanmean(mean_age,3));

%--- Summer 2005 ---%
stind=find(time==datenum(2005,6,1))-1;
enind=find(time==datenum(2005,9,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_06'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
mean_age=squeeze(double(var2./var1));
mean_age(find(var1<1e-6))=nan;

summ.ma(1,:,:)=squeeze(nanmean(mean_age,3));

clear var1 var2 mean_age;

nanmask=nanmaskh;

%--- Autumn 2005 ---%
stind=find(time==datenum(2005,9 ,1))-1;
enind=find(time==datenum(2005,12,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_06'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
mean_age=squeeze(double(var2./var1));
mean_age(find(var1<1e-6))=nan;

autumn.ma(1,:,:)=squeeze(nanmean(mean_age,3));

clear var1 var2 mean_age;

nanmask=nanmaskh;

%--- Winter 2005/6 ---%
stind=find(time==datenum(2005,12,1))-1;
enind=find(time==datenum(2006, 3,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_06'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
mean_age=squeeze(double(var2./var1));
mean_age(find(var1<1e-6))=nan;

winter.ma(1,:,:)=squeeze(nanmean(mean_age,3));

clear var1 var2 mean_age;

nanmask=nanmaskh;

%--- Spring 2006 ---%
stind=find(time==datenum(2006,3,1))-1;
enind=find(time==datenum(2006,6,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_06'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
mean_age=squeeze(double(var2./var1));
mean_age(find(var1<1e-6))=nan;

spring.ma(2,:,:)=squeeze(nanmean(mean_age,3));

clear var1 var2 mean_age;

nanmask=nanmaskh;

%--- Summer 2006 ---%
stind=find(time==datenum(2006,6,1))-1;
enind=find(time==datenum(2006,9,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_06'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
mean_age=squeeze(double(var2./var1));
mean_age(find(var1<1e-6))=nan;

summ.ma(2,:,:)=squeeze(nanmean(mean_age,3));

clear var1 var2 mean_age;

nanmask=nanmaskh;

%--- Autumn 2006 ---%
stind=find(time==datenum(2006,9 ,1))-1;
enind=find(time==datenum(2006,12,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_06'),...
                  [0 0 35 stind],[402 482 1 length(stind:enind-1)]);
mean_age=squeeze(double(var2./var1));
mean_age(find(var1<1e-6))=nan;

autumn.ma(2,:,:)=squeeze(nanmean(mean_age,3));

clear var1 var2 mean_age;

nanmask=nanmaskh;

%========================================================================%
%=                       Plot Spring 2005 Mean Age                      =%
%========================================================================%
figure(1); clf;
ax1=axes('position',[0.01 0.05 0.8 0.93]);
m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                 'lat' ,[min(min(lat)) max(max(lat))]);
m_grid('tickdir','on');
hold on;
m_pcolor(lon,lat,squeeze(spring.ma(1,:,:))'.*nanmaskh); shading flat;
m_contour(lon,lat,hmask',[0 0],'k');
m_text(-95.5,38,'Mean water age','fontsize',14,'fontweight','bold');
m_text(-92.5,36.75,'(surface)','fontsize',14,'fontweight','bold');
caxis([0 365]);

%--- Add Colorbar (scaled by current day) ---%
cax=axes('position',[0.82 0.075 0.075 0.7]); axis off;
caxis([0 365]);
hCbar     =colorbar('position',[0.82 0.075 0.075 0.7]);
hCbarT=text(-0.13,-0.05,'Age (days)','fontsize',12,...
                                     'fontname','helvetica',...
                                     'fontweight','bold');
set(hCbar,'linewidth',1,'ytick',[0:25:(floor(365/25)*25)-1 365],...
          'fontsize',10,'fontname','helvetica','fontweight','bold');


%--- Make set of text at top for model run info ---%
tax=axes('position',[0.7 0.8 0.275 0.175]); axis off;
hmodelname=text(-0.075,0.95,'US East Coast water age model',...
                           'fontweight','bold','fontsize',12);
hNCSU     =text(0.3,0.775,'NC State Univeristy (2014)','fontsize',10);             
hmodelres =text(0.645,0.56,'\Delta x,\Delta y = 1/10^{o}','fontsize',10);          
hmodelvert=text(0.865,0.395,'k = 36','fontsize',10);
hmodeldate=text(1.06,0.2,'Mar - May 2005',...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');

set(gcf,'color','w');
wysiwyg;
export_fig /home/actodd/PLOTS/USeast/age/spring2005age.png -png -r150 -painters


%========================================================================%
%=                       Plot Summer 2005 Mean Age                      =%
%========================================================================%
figure(1); clf;
ax1=axes('position',[0.01 0.05 0.8 0.93]);
m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                 'lat' ,[min(min(lat)) max(max(lat))]);
m_grid('tickdir','on');
hold on;
m_pcolor(lon,lat,squeeze(summ.ma(1,:,:))'.*nanmaskh); shading flat;
m_contour(lon,lat,hmask',[0 0],'k');
m_text(-95.5,38,'Mean water age','fontsize',14,'fontweight','bold');
m_text(-92.5,36.75,'(surface)','fontsize',14,'fontweight','bold');
caxis([0 365]);

%--- Add Colorbar (scaled by current day) ---%
cax=axes('position',[0.82 0.075 0.075 0.7]); axis off;
caxis([0 365]);
hCbar     =colorbar('position',[0.82 0.075 0.075 0.7]);
hCbarT=text(-0.13,-0.05,'Age (days)','fontsize',12,...
                                     'fontname','helvetica',...
                                     'fontweight','bold');
set(hCbar,'linewidth',1,'ytick',[0:25:(floor(365/25)*25)-1 365],...
          'fontsize',10,'fontname','helvetica','fontweight','bold');


%--- Make set of text at top for model run info ---%
tax=axes('position',[0.7 0.8 0.275 0.175]); axis off;
hmodelname=text(-0.075,0.95,'US East Coast water age model',...
                           'fontweight','bold','fontsize',12);
hNCSU     =text(0.3,0.775,'NC State Univeristy (2014)','fontsize',10);             
hmodelres =text(0.645,0.56,'\Delta x,\Delta y = 1/10^{o}','fontsize',10);          
hmodelvert=text(0.865,0.395,'k = 36','fontsize',10);
hmodeldate=text(1.06,0.2,'Jun - Aug 2005',...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');

set(gcf,'color','w');
wysiwyg;
export_fig /home/actodd/PLOTS/USeast/age/summer2005age.png -png -r150 -painters


%========================================================================%
%=                       Plot Autumn 2005 Mean Age                      =%
%========================================================================%
figure(1); clf;
ax1=axes('position',[0.01 0.05 0.8 0.93]);
m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                 'lat' ,[min(min(lat)) max(max(lat))]);
m_grid('tickdir','on');
hold on;
m_pcolor(lon,lat,squeeze(autumn.ma(1,:,:))'.*nanmaskh); shading flat;
m_contour(lon,lat,hmask',[0 0],'k');
m_text(-95.5,38,'Mean water age','fontsize',14,'fontweight','bold');
m_text(-92.5,36.75,'(surface)','fontsize',14,'fontweight','bold');
caxis([0 365]);

%--- Add Colorbar (scaled by current day) ---%
cax=axes('position',[0.82 0.075 0.075 0.7]); axis off;
caxis([0 365]);
hCbar     =colorbar('position',[0.82 0.075 0.075 0.7]);
hCbarT=text(-0.13,-0.05,'Age (days)','fontsize',12,...
                                     'fontname','helvetica',...
                                     'fontweight','bold');
set(hCbar,'linewidth',1,'ytick',[0:25:(floor(365/25)*25)-1 365],...
          'fontsize',10,'fontname','helvetica','fontweight','bold');


%--- Make set of text at top for model run info ---%
tax=axes('position',[0.7 0.8 0.275 0.175]); axis off;
hmodelname=text(-0.075,0.95,'US East Coast water age model',...
                           'fontweight','bold','fontsize',12);
hNCSU     =text(0.3,0.775,'NC State Univeristy (2014)','fontsize',10);             
hmodelres =text(0.645,0.56,'\Delta x,\Delta y = 1/10^{o}','fontsize',10);          
hmodelvert=text(0.865,0.395,'k = 36','fontsize',10);
hmodeldate=text(1.06,0.2,'Sep - Nov 2005',...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');

set(gcf,'color','w');
wysiwyg;
export_fig /home/actodd/PLOTS/USeast/age/autumn2005age.png -png -r150 -painters


%========================================================================%
%=                       Plot Winter 2005 Mean Age                      =%
%========================================================================%
figure(1); clf;
ax1=axes('position',[0.01 0.05 0.8 0.93]);
m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                 'lat' ,[min(min(lat)) max(max(lat))]);
m_grid('tickdir','on');
hold on;
m_pcolor(lon,lat,squeeze(winter.ma(1,:,:))'.*nanmaskh); shading flat;
m_contour(lon,lat,hmask',[0 0],'k');
m_text(-95.5,38,'Mean water age','fontsize',14,'fontweight','bold');
m_text(-92.5,36.75,'(surface)','fontsize',14,'fontweight','bold');
caxis([0 365]);

%--- Add Colorbar (scaled by current day) ---%
cax=axes('position',[0.82 0.075 0.075 0.7]); axis off;
caxis([0 365]);
hCbar     =colorbar('position',[0.82 0.075 0.075 0.7]);
hCbarT=text(-0.13,-0.05,'Age (days)','fontsize',12,...
                                     'fontname','helvetica',...
                                     'fontweight','bold');
set(hCbar,'linewidth',1,'ytick',[0:25:(floor(365/25)*25)-1 365],...
          'fontsize',10,'fontname','helvetica','fontweight','bold');


%--- Make set of text at top for model run info ---%
tax=axes('position',[0.7 0.8 0.275 0.175]); axis off;
hmodelname=text(-0.075,0.95,'US East Coast water age model',...
                           'fontweight','bold','fontsize',12);
hNCSU     =text(0.3,0.775,'NC State Univeristy (2014)','fontsize',10);             
hmodelres =text(0.645,0.56,'\Delta x,\Delta y = 1/10^{o}','fontsize',10);          
hmodelvert=text(0.865,0.395,'k = 36','fontsize',10);
hmodeldate=text(1.06,0.2,'Dec 2005 - Feb 2006',...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');

set(gcf,'color','w');
wysiwyg;
export_fig /home/actodd/PLOTS/USeast/age/winter2005age.png -png -r150 -painters


%========================================================================%
%=                       Plot Spring 2006 Mean Age                      =%
%========================================================================%
figure(1); clf;
ax1=axes('position',[0.01 0.05 0.8 0.93]);
m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                 'lat' ,[min(min(lat)) max(max(lat))]);
m_grid('tickdir','on');
hold on;
m_pcolor(lon,lat,squeeze(spring.ma(2,:,:))'.*nanmaskh); shading flat;
m_contour(lon,lat,hmask',[0 0],'k');
m_text(-95.5,38,'Mean water age','fontsize',14,'fontweight','bold');
m_text(-92.5,36.75,'(surface)','fontsize',14,'fontweight','bold');
caxis([0 365]);

%--- Add Colorbar (scaled by current day) ---%
cax=axes('position',[0.82 0.075 0.075 0.7]); axis off;
caxis([0 365]);
hCbar     =colorbar('position',[0.82 0.075 0.075 0.7]);
hCbarT=text(-0.13,-0.05,'Age (days)','fontsize',12,...
                                     'fontname','helvetica',...
                                     'fontweight','bold');
set(hCbar,'linewidth',1,'ytick',[0:25:(floor(365/25)*25)-1 365],...
          'fontsize',10,'fontname','helvetica','fontweight','bold');


%--- Make set of text at top for model run info ---%
tax=axes('position',[0.7 0.8 0.275 0.175]); axis off;
hmodelname=text(-0.075,0.95,'US East Coast water age model',...
                           'fontweight','bold','fontsize',12);
hNCSU     =text(0.3,0.775,'NC State Univeristy (2014)','fontsize',10);             
hmodelres =text(0.645,0.56,'\Delta x,\Delta y = 1/10^{o}','fontsize',10);          
hmodelvert=text(0.865,0.395,'k = 36','fontsize',10);
hmodeldate=text(1.06,0.2,'Mar - May 2006',...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');

set(gcf,'color','w');
wysiwyg;
export_fig /home/actodd/PLOTS/USeast/age/spring2006age.png -png -r150 -painters


%========================================================================%
%=                       Plot Summer 2006 Mean Age                      =%
%========================================================================%
figure(1); clf;
ax1=axes('position',[0.01 0.05 0.8 0.93]);
m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                 'lat' ,[min(min(lat)) max(max(lat))]);
m_grid('tickdir','on');
hold on;
m_pcolor(lon,lat,squeeze(summ.ma(2,:,:))'.*nanmaskh); shading flat;
m_contour(lon,lat,hmask',[0 0],'k');
m_text(-95.5,38,'Mean water age','fontsize',14,'fontweight','bold');
m_text(-92.5,36.75,'(surface)','fontsize',14,'fontweight','bold');
caxis([0 365]);

%--- Add Colorbar (scaled by current day) ---%
cax=axes('position',[0.82 0.075 0.075 0.7]); axis off;
caxis([0 365]);
hCbar     =colorbar('position',[0.82 0.075 0.075 0.7]);
hCbarT=text(-0.13,-0.05,'Age (days)','fontsize',12,...
                                     'fontname','helvetica',...
                                     'fontweight','bold');
set(hCbar,'linewidth',1,'ytick',[0:25:(floor(365/25)*25)-1 365],...
          'fontsize',10,'fontname','helvetica','fontweight','bold');


%--- Make set of text at top for model run info ---%
tax=axes('position',[0.7 0.8 0.275 0.175]); axis off;
hmodelname=text(-0.075,0.95,'US East Coast water age model',...
                           'fontweight','bold','fontsize',12);
hNCSU     =text(0.3,0.775,'NC State Univeristy (2014)','fontsize',10);             
hmodelres =text(0.645,0.56,'\Delta x,\Delta y = 1/10^{o}','fontsize',10);          
hmodelvert=text(0.865,0.395,'k = 36','fontsize',10);
hmodeldate=text(1.06,0.2,'Jun - Aug 2006',...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');

set(gcf,'color','w');
wysiwyg;
export_fig /home/actodd/PLOTS/USeast/age/summer2006age.png -png -r150 -painters


%========================================================================%
%=                       Plot Autumn 2006 Mean Age                      =%
%========================================================================%
figure(1); clf;
ax1=axes('position',[0.01 0.05 0.8 0.93]);
m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                 'lat' ,[min(min(lat)) max(max(lat))]);
m_grid('tickdir','on');
hold on;
m_pcolor(lon,lat,squeeze(autumn.ma(2,:,:))'.*nanmaskh); shading flat;
m_contour(lon,lat,hmask',[0 0],'k');
m_text(-95.5,38,'Mean water age','fontsize',14,'fontweight','bold');
m_text(-92.5,36.75,'(surface)','fontsize',14,'fontweight','bold');
caxis([0 365]);

%--- Add Colorbar (scaled by current day) ---%
cax=axes('position',[0.82 0.075 0.075 0.7]); axis off;
caxis([0 365]);
hCbar     =colorbar('position',[0.82 0.075 0.075 0.7]);
hCbarT=text(-0.13,-0.05,'Age (days)','fontsize',12,...
                                     'fontname','helvetica',...
                                     'fontweight','bold');
set(hCbar,'linewidth',1,'ytick',[0:25:(floor(365/25)*25)-1 365],...
          'fontsize',10,'fontname','helvetica','fontweight','bold');


%--- Make set of text at top for model run info ---%
tax=axes('position',[0.7 0.8 0.275 0.175]); axis off;
hmodelname=text(-0.075,0.95,'US East Coast water age model',...
                           'fontweight','bold','fontsize',12);
hNCSU     =text(0.3,0.775,'NC State Univeristy (2014)','fontsize',10);             
hmodelres =text(0.645,0.56,'\Delta x,\Delta y = 1/10^{o}','fontsize',10);          
hmodelvert=text(0.865,0.395,'k = 36','fontsize',10);
hmodeldate=text(1.06,0.2,'Sep - Nov 2006',...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');

set(gcf,'color','w');
wysiwyg;
export_fig /home/actodd/PLOTS/USeast/age/autumn2006age.png -png -r150 -painters

