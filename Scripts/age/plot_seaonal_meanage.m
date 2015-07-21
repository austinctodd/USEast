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
addpath(genpath('/Users/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=============================================================================%
%=                                                                           =%
%=       !!!!!!!!!!!        USER-DEFINED VARIABLES      !!!!!!!!!!           =%
%=                                                                           =%
%=============================================================================%
tracer='shelf'; % 'shelf' or 'deep' or 'river'
pdir='/Volumes/Black_box/Data/USeast-age/output/clim/averages/';
plotdir='/Users/actodd/Documents/Work/Projects/US East/FIGS/';

%=============================================================================%
%=                  Setup some grid and preliminary arrays                   =%
%=============================================================================%
refdate=datenum(2004,01,01); % Why the hell do we choose this date???

gridfile  = '/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc'

ncid=netcdf.open(gridfile,'nowrite');
  h    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'         ));
  lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'   ));
  lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'   ));
  mask1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_shelf'));
  hmask=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'  ));
netcdf.close(ncid)

lon=lon'; lat=lat';
nanmask1=mask1'./mask1';
nanmaskh=hmask'./hmask';
time=time./86400+refdate;

%=============================================================================%
%=                   Load seasonal age data from MAT file                    =%
%=============================================================================%
input_file  = [pdir,'avg_3hrly.nc'];

ncid=netcdf.open(input_file,'nowrite');
  age_01=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mean_age_01'),[],[]);
  age_02=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mean_age_02'));
  age_03=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mean_age_03'));
  time =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'));
netcdf.close(ncid)


%==============================================================================%
%=                                Plot Mean Age                               =%
%==============================================================================%
for i=1:4
  if i==1
    var=squeeze(spring.ma(1,:,:))';
    monthlab='Mar - May 2005';
    outfname='spring2005age';
  elseif i==2
    var=squeeze(summ.ma(1,:,:))';
    monthlab='Jun - Aug 2005';
    outfname='summer2005age';
  elseif i==3
    var=squeeze(autumn.ma(1,:,:))';
    monthlab='Sep - Nov 2005';
    outfname='autumn2005age';
  elseif i==4
    var=squeeze(winter.ma(1,:,:))';
    monthlab='Dec 2005 - Feb 2006';
    outfname='winter2005age';
  elseif i==5
    var=squeeze(spring.ma(2,:,:))';
    monthlab='Mar - May 2006';
    outfname='spring2006age';
  elseif i==6
    var=squeeze(summ.ma(2,:,:))';
    monthlab='Jun - Aug 2006';
    outfname='summer2006age';
  else
    var=squeeze(autumn.ma(2,:,:))';
    monthlab='Sep - Nov 2006';
    outfname='autumn2006age';
  end

  figure(1); clf;
  ax1=axes('position',[0.01 0.05 0.8 0.93]);
  m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                   'lat' ,[min(min(lat)) max(max(lat))]);
  m_grid('tickdir','on');
  hold on;
  m_pcolor(lon,lat,var.*nanmaskh); shading flat;
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
  hmodeldate=text(1.06,0.2,monthlab,...
                      'fontsize',14,'fontweight','bold',...
                      'horizontalalignment','right');

  set(gcf,'color','w');
  wysiwyg;
  print('-depsc2','-painters',[plotdir,outfname,'.eps']);
  %export_fig /home/actodd/PLOTS/USeast/age/autumn2006age.png -png -r150 -painters

end
  