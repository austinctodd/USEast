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
plotdir='/home/actodd/PLOTS/USeast/age/';

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
netcdf.close(ncid)

lon=lon'; lat=lat';
nanmask1=mask1'./mask1';
nanmaskh=hmask'./hmask';
time=time./86400+refdate;
%=============================================================================%
%=                   Load seasonal age data from MAT file                    =%
%=============================================================================%
load meanrt.mat

%==============================================================================%
%=                                Plot Mean Age                               =%
%==============================================================================%
for i=1:7
  if i==1
    var=squeeze(spring.rt(1,:,:))';
    monthlab='Mar - May 2005';
    outfname='spring2005rt';
  elseif i==2
    var=squeeze(summ.rt(1,:,:))';
    monthlab='Jun - Aug 2005';
    outfname='summer2005rt';
  elseif i==3
    var=squeeze(autumn.rt(1,:,:))';
    monthlab='Sep - Nov 2005';
    outfname='autumn2005rt';
  elseif i==4
    var=squeeze(winter.rt(1,:,:))';
    monthlab='Dec 2005 - Feb 2006';
    outfname='winter2005rt';
  elseif i==5
    var=squeeze(spring.rt(2,:,:))';
    monthlab='Mar - May 2006';
    outfname='spring2006rt';
  elseif i==6
    var=squeeze(summ.rt(2,:,:))';
    monthlab='Jun - Aug 2006';
    outfname='summer2006rt';
  else
    var=squeeze(autumn.rt(2,:,:))';
    monthlab='Sep - Nov 2006';
    outfname='autumn2006rt';
  end

  figure(1); clf;
  ax1=axes('position',[0.01 0.05 0.8 0.93]);
  m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                   'lat' ,[min(min(lat)) max(max(lat))]);
  m_grid('tickdir','on');
  hold on;
  m_pcolor(lon,lat,var.*nanmask1); shading flat;
  m_contour(lon,lat,hmask',[0 0],'k');
  %m_text(-95.5,38,'Mean residence','fontsize',14,'fontweight','bold');
  %m_text(-94.8,36.75,'time (surface)','fontsize',14,'fontweight','bold');
  caxis([0 250]);

  %--- Add Colorbar (scaled by current day) ---%
  cax=axes('position',[0.82 0.075 0.075 0.7]); axis off;
  day=250;
  caxis([0 day]);
  hCbar     =colorbar('position',[0.82 0.075 0.075 0.7]);
  hCbarT=text(1.83,0.73,'Residence time (days)','fontsize',12,...
                                     'fontname','helvetica',...
                                     'fontweight','bold','Rotation',-90);
  set(hCbar,'linewidth',1,'ytick',[0:25:(floor(day/25)*25)-1 day],...
            'fontsize',10,'fontname','helvetica','fontweight','bold');


  %--- Make set of text at top for model run info ---%
  tax=axes('position',[0.7 0.8 0.275 0.175]); axis off;
  hmodelname=text(-0.075,0.95,'US East Coast water age model',...
                              'fontweight','bold','fontsize',12);
  hNCSU     =text(0.3,0.775,'NC State Univeristy (2014)','fontsize',10);             
  hmodelres =text(0.645,0.56,'\Delta x,\Delta y = 1/10^{o}','fontsize',10);          
  hmodelvert=text(0.865,0.395,'k = 36','fontsize',10);
  hmodeldate=text(1.06,0.2,monthlab,...
                      'fontsize',14,'fontweight','bold',...
                      'horizontalalignment','right');

  set(gcf,'color','w');
  wysiwyg;
  print('-depsc2','-painters',[plotdir,outfname,'.eps']);
  eval(['export_fig ',[plotdir,outfname],'.png -png -r150 -painters']);

end
  