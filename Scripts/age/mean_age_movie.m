%=============================================================================%
%=                                                                           =%
%=  Created by Austin C. Todd, NCSU, 2013                                    =%
%=                                                                           =%
%=  Program is designed to read in age information from the US East model    =%
%=  and calculate the mean age.  The model has three different "source"      =%
%=  water regions:  1. Shelf water, 2. River water, 3. Deep ocean water.     =%
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
addpath('/home/actodd/MYMATLAB/');
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));
addpath(genpath('/home/actodd/MYMATLAB/m_map/'));
addpath(genpath('/he_data/he/zhigang/Matlab_Codes/mexcdf2009/'))
format long g; format compact;

clear all;
%=============================================================================%
%=                                                                           =%
%=       !!!!!!!!!!!        USER-DEFINED VARIABLES      !!!!!!!!!!           =%
%=                                                                           =%
%=============================================================================%
tracer='deep'; % 'shelf' or 'deep' or 'river'
pdir='/gpfs_share/rhe/actodd/USeast-age/output/';
%pdir='/gpfs_share/rhe/actodd/USeast-age/output';

%=============================================================================%
%=                  Setup some grid and preliminary arrays                   =%
%=============================================================================%
refdate=datenum(2004,1,1); 

gridfile=['/he_data/he/actodd/ROMS-age/Data/frc/USeast-grid.nc'];
h      =nc_varget(gridfile,'h'         );
lon    =nc_varget(gridfile,'lon_rho'   );
lat    =nc_varget(gridfile,'lat_rho'   );
smask  =nc_varget(gridfile,'mask_shelf' );
hmask  =nc_varget(gridfile,'mask_rho'  );
snanmask=smask./smask;
hnanmask =hmask./hmask;

%=============================================================================%
%=        Load in colorbar and make one that is just shades of purple        =%
%=============================================================================%
%cmap=load('/raid0/actodd/LTRANS-utils/density_colormap.txt');
%cmap(end,:)=cmap(end-1,:);
%for i=29:-1:1
%  cmap(i,1)=cmap(i+1,1)+(0.9-cmap(30,1))/30;
%  cmap(i,2)=cmap(i+1,2)+(0.9/30);
%  cmap(i,3)=cmap(i+1,3)+(0.9-cmap(30,3))/30;
%end
%cmap=cmap(end:-1:1,:);
%cmap(end,:)=[1 1 1];

%=============================================================================%
%=           Loop to read in files and plot surface age data                 =%
%=============================================================================%
count=1;
for day=366:366

  if (day < 93)
    pdir='/gpfs_share/rhe/actodd/USeast-age/output/OBC/';
  else
    pdir='/gpfs_share/rhe/actodd/USeast-age/output/';
  end
  if day==1
    filenum=1; ind=0;
  elseif day==2
    filenum=1; ind=8;
  else
    filenum=day-1; ind=7;
  end
  fname=[pdir,'useast_his_',sprintf('%04i',filenum),'.nc'];

  disp(['Day : ',sprintf('%02i',day)]);

  if strcmp(tracer,'shelf')
    var1=nc_varget(fname,['age_01'],[ind 35 0 0],[1 1 482 402]);
    var2=nc_varget(fname,['age_02'],[ind 35 0 0],[1 1 482 402]);
    nanmask=nanmask1;
  elseif strcmp(tracer,'deep');
    var1=nc_varget(fname,['age_05'],[ind 35 0 0],[1 1 482 402]);
    var2=nc_varget(fname,['age_06'],[ind 35 0 0],[1 1 482 402]);
    hnanmask=hnanmask;
  else
    var1=nc_varget(fname,['age_03'],[ind 35 0 0],[1 1 482 402]);
    var2=nc_varget(fname,['age_04'],[ind 35 0 0],[1 1 482 402]);
    nanmask=nanmask2;
  end
  %var1(find(var1 < 1e-5))=0; 
  %var2(find(var1 < 1e-5))=0;
  mean_age=var2./var1;
  mean_age(find(var1<1e-6))=nan;
  %mean_age(find(mean_age < 0 ))=nan;
  %mean_age(find(mean_age >60))=60;  
 
  %--- PLOT Age on Map---%
  clf;
  ax1=axes('position',[0.01 0.05 0.8 0.93]);
  m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
              'lat' ,[min(min(lat)) max(max(lat))]);
  m_grid('tickdir','on');
  hold on; 
  m_pcolor(lon,lat,mean_age.*hnanmask); shading flat;
  m_contour(lon,lat,hmask,[0 0],'k');
  m_text(-95.5,38,'Mean water age','fontsize',14,'fontweight','bold');
  m_text(-92.5,36.75,'(surface)','fontsize',14,'fontweight','bold');
  caxis([0 day]);
  
  %--- Add Colorbar (scaled by current day) ---%
  cax=axes('position',[0.82 0.075 0.075 0.7]); axis off;
  cbarlength=day/365; caxis([0 day]);
  hCbar =colorbar('position',[0.82 0.075 0.075 0.7*cbarlength]);
  hCbarT=text(-0.13,-0.05,'Age (days)','fontsize',12,...
                                       'fontname','helvetica',...
                                       'fontweight','bold');
  set(hCbar,'linewidth',1,'ytick',[0:25:(floor(day/25)*25)-1 day],...
            'fontsize',10,'fontname','helvetica','fontweight','bold');
  
                                   
  %--- Make set of text at top for model run info ---%
  tax=axes('position',[0.7 0.8 0.275 0.175]); axis off;
  hmodelname=text(-0.075,0.95,'US East Coast water age model',...
                             'fontweight','bold','fontsize',12);
  hNCSU     =text(0.3,0.775,'NC State Univeristy (2014)','fontsize',10);                         
  hmodelres =text(0.645,0.56,'\Delta x,\Delta y = 1/12^{o}','fontsize',10);                       
  hmodelvert=text(0.865,0.395,'k = 36','fontsize',10);                       
  hmodeldate=text(1.06,0.2,datestr(datenum(2004,01,day),1),...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');
  
  set(gcf,'color','w');
  
  drawnow;
  wysiwyg;
  pause(0.25); 
  eval(['export_fig /home/actodd/PLOTS/USeast/age/age_',...
        sprintf('%03i',day),'.png -png -r150 -painters']);
%  print('-depsc2','-painters',['/raid0/actodd/PLOTS/USeast-age/frames/',...
%                               'meanage_',sprintf('%02i',day),'.eps']);

end


