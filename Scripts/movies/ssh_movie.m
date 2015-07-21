%=============================================================================%
%=                                                                           =%
%=  Created by Austin C. Todd, NCSU, 2013                                    =%
%=                                                                           =%
%=  Program is designed to read in SSH from the US East model and plot    =%
%=  and calculate the mean age.  The model has three different "source"      =%
%=  water regions:  1. Shelf water, 2. River water, 3. Deep ocean water.     =%
%=                                                                           =%
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
for day=1:8:366*8

  pdir='/gpfs_share/rhe/actodd/USeast-age/output/';
  fname=[pdir,'useast_his.nc'];

  disp(['Day : ',sprintf('%02i',(day-1)/8+1)]);
  
  ssh=nc_varget(fname,['zeta'],[day-1 0 0],[1 482 402]);
  nanmask=hnanmask;
 
  %--- PLOT SSH on Map---%
  clf;
  ax1=axes('position',[0.01 0.05 0.8 0.93]);
  m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
              'lat' ,[min(min(lat)) max(max(lat))]);
  m_grid('tickdir','on');
  hold on; 
  m_pcolor(lon,lat,ssh.*hnanmask); shading flat;
  m_contour(lon,lat,hmask,[0 0],'k');
  m_text(-92.5,36.75,'Anomoly (m)','fontsize',14,'fontweight','bold');
  m_text(-95.5,38,'Sea Surface Height','fontsize',14,'fontweight','bold');
  caxis([-1.5 1.5]);
  
  %--- Add Colorbar (scaled by current day) ---%
  cax=axes('position',[0.82 0.075 0.075 0.7]); axis off;
  cbarlength=1; caxis([-1.5 1.5]);
  hCbar =colorbar('position',[0.82 0.075 0.075 0.7]);
  hCbarT=text(0,-0.04,'SSH (m)','fontsize',12,...
                                       'fontname','helvetica',...
                                       'fontweight','bold');
  set(hCbar,'linewidth',1,'ytick',[-1.5:0.25:1.5],...
            'fontsize',10,'fontname','helvetica','fontweight','bold');
  
                                   
  %--- Make set of text at top for model run info ---%
  tax=axes('position',[0.7 0.8 0.275 0.175]); axis off;
  hmodelname=text(-0.075,0.95,'US East Coast water age model',...
                             'fontweight','bold','fontsize',12);
  hNCSU     =text(0.3,0.775,'NC State Univeristy (2014)','fontsize',10);                         
  hmodelres =text(0.645,0.56,'\Delta x,\Delta y = 1/12^{o}','fontsize',10);                       
  hmodelvert=text(0.865,0.395,'k = 36','fontsize',10);                       
  hmodeldate=text(1.06,0.2,datestr(datenum(2005,01,(day-1)./8+1),1),...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');
  
  set(gcf,'color','w');
  
  wysiwyg;
  eval(['export_fig /home/actodd/PLOTS/USeast/frames/ssh_',...
        sprintf('%03i',day),'.png -png -r150 -painters']);

end


