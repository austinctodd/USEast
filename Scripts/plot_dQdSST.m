%============================= plot_eke.m ================================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program is designed to read in the KE and EKE from .mat files & =%
%=  plot maps of mean KE and mean EKE. For the KE calculation, the       =%
%=  velocities are broken down into u=u_ + u', where u_ is the time mean =%
%=  velocity, and u' is the time-varying perturbation velocity. The KE   =%
%=  is then KE=0.5*(u_^2 + v_^2) and EKE=0.5*(u'^2 + v'^2).              =%
%=                                                                       =%
%=========================================================================%

%--- Add various libraries and paths ---%
addpath(genpath('/Users/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/Users/actodd/MYMATLAB/');
addpath(genpath('/Users/actodd/MYMATLAB/m_map'));
format long g; format compact;
set(0,'defaultfigurerenderer','zbuffer');

%=========================================================================%
%=                  Set file input/output directories                    =%
%=========================================================================%
plotdir='/Users/actodd/PLOTS/USeast/temp/frames/';
datadir='/Volumes/Black_box/Data/USeast/Data/atmos/';
climname='/Volumes/Black_box/Data/USeast/Data/clim/USeast-clim-monthly.nc';

%=========================================================================%
%=                      Load ROMS grid information                       =%
%=========================================================================%
disp('Reading in ROMS grid information');

fname=['/Volumes/Black_box/Data/USeast-age/output/clim/testdQdSST/',...
       'ERA/useast_his.nc'];

ncid=netcdf.open(fname,'nowrite');
 lon        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'    ),...
                                                [0 0],[402-25 482]);
 lat        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'    ),...
                                                [0 0],[402-25 482]);
 mask       =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'   ),...
                                                [0 0],[402-25 482]);
 h          =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'          ),...
                                                [0 0],[402-25 482]);
netcdf.close(ncid);

[sizelt,sizeln]=size(mask);
nanmask=mask./mask;

%=========================================================================%
%=                         Load ERA temps/dQdSST                         =%
%=========================================================================%
ncid=netcdf.open([datadir,'ERA-dQdSST.nc'],'nowrite');
  era.lat   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
  era.lon   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
  era.dQdSST_raw=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dQdSST'));
  era.sst_raw   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'SST'));
netcdf.close(ncid);
era.lat   =double(era.lat);
era.lon   =double(era.lon);
era.dQdSST_raw=double(squeeze(era.dQdSST_raw));
era.sst_raw   =double(squeeze(era.sst_raw));

[x,y]=meshgrid(era.lon,era.lat);

%--- Interpolate to ROMS grid ---%
for i=1:12
    era.dQdSST(:,:,i)=griddata(x,y,squeeze(era.dQdSST_raw(:,:,i))',...
                               lon,lat);
    era.sst(   :,:,i)=griddata(x,y,squeeze(era.sst_raw(:,:,i))',...
                               lon,lat);
end

%=========================================================================%
%=                        Load CFSR temps/dQdSST                         =%
%=========================================================================%
ncid=netcdf.open([datadir,'CFSR-dQdSST.nc'],'nowrite');
  cfsr.lat   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'),...
                                          [0 0],[402-25 482]);
  cfsr.lon   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'),...
                                          [0 0],[402-25 482]);
  cfsr.dQdSST=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dQdSST'),...
                                          [0 0 0],[402-25 482 12]);
  cfsr.sst=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'SST'),...
                                          [0 0 0],[402-25 482 12]);
netcdf.close(ncid);
cfsr.lat   =double(cfsr.lat);
cfsr.lon   =double(cfsr.lon);
cfsr.dQdSST=double(squeeze(cfsr.dQdSST));
cfsr.sst   =double(squeeze(cfsr.sst   ));

%=========================================================================%
%=                      Load ERA/HYCOM temps/dQdSST                      =%
%=========================================================================%
ncid=netcdf.open([datadir,'ERA-dQdSST-hycom.nc'],'nowrite');
  era_hycom.lat   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'),...
                                          [0 0],[402-25 482]);
  era_hycom.lon   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'),...
                                          [0 0],[402-25 482]);
  era_hycom.dQdSST=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dQdSST'),...
                                          [0 0 0],[402-25 482 12]);
  era_hycom.sst=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'SST'),...
                                          [0 0 0],[402-25 482 12]);
netcdf.close(ncid);
era_hycom.lat   =double(era_hycom.lat);
era_hycom.lon   =double(era_hycom.lon);
era_hycom.dQdSST=double(squeeze(era_hycom.dQdSST));
era_hycom.sst   =double(squeeze(era_hycom.sst   ));

%=========================================================================%
%=                     Load CFSR/HYCOM temps/dQdSST                      =%
%=========================================================================%
ncid=netcdf.open([datadir,'CFSR-dQdSST-hycom.nc'],'nowrite');
  cfsr_hycom.lat   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'),...
                                          [0 0],[402-25 482]);
  cfsr_hycom.lon   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'),...
                                          [0 0],[402-25 482]);
  cfsr_hycom.dQdSST=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dQdSST'),...
                                          [0 0 0],[402-25 482 12]);
  cfsr_hycom.sst=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'SST'),...
                                          [0 0 0],[402-25 482 12]);
netcdf.close(ncid);
cfsr_hycom.lat   =double(cfsr_hycom.lat);
cfsr_hycom.lon   =double(cfsr_hycom.lon);
cfsr_hycom.dQdSST=double(squeeze(cfsr_hycom.dQdSST));
cfsr_hycom.sst   =double(squeeze(cfsr_hycom.sst   ));

%=========================================================================%
%=                     Set different axes properties                     =%
%=========================================================================%
ax1pos(1,:)=[0.035 0.53 0.35 0.465];
ax1pos(2,:)=[0.340 0.53 0.35 0.465];
ax1pos(3,:)=[0.035 0.03 0.35 0.465];
ax1pos(4,:)=[0.340 0.03 0.35 0.465];
ax1pos(5,:)=[0.645 0.53 0.35 0.465];

exnames={'era';'era_hycom';'cfsr';'cfsr_hycom'};
extitle={'ERA';'ERA/HYCOM';'CFSR';'CFSR/HYCOM'};

figure(1); clf;
return;
for ex=1:4
    
  eval(['temp=squeeze(nanmean(',exnames{ex},'.dQdSST(:,:,10:12),3));']);
  
  %=======================================================================%
  %=                   Main plot panel for each exp                      =%
  %=======================================================================%
  eval(['ax',sprintf('%1i',ex),'=axes(''position'',ax1pos(',...
             sprintf('%1i',ex),',:));']);
  m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                    'lat',[min(min(lat)) max(max(lat))]);
  if (ex==1)
    m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
           'xticklabels',char({'96^{o}W';'88^{o}W';'80^{o}W';'72^{o}W';''}),...
	       'yticklabels',char({'16^{o}N ';'24^{o}N ';'32^{o}N ';'40^{o}N ';}),...
	       'fontsize',10)
    elseif (ex==2)
      m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
             'xticklabels',char({'96^{o}W';'88^{o}W';'80^{o}W';'72^{o}W';''}),...
	     'yticklabels',char({'';'';'32^{o}N ';'40^{o}N '}),...
	     'fontsize',10)
    elseif (ex==3)
      m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
             'xticklabels',char({'96^{o}W';'88^{o}W';'80^{o}W';'72^{o}W';''}),...
	     'yticklabels',char({'16^{o}N ';'24^{o}N ';'32^{o}N ';'40^{o}N ';}),...
	     'fontsize',10)
    else
      m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
             'xticklabels',char({'96^{o}W';'88^{o}W';'80^{o}W';'72^{o}W';'64^{o}W'}),...
	     'yticklabels',char({'';'';'32^{o}N ';'40^{o}N ';}),...
	     'fontsize',10)
    end
    hold on;
    m_contour(lon,lat,mask,[0 0],'k','linewidth',1.5);
    pcc=m_pcolor(lon,lat,-temp.*mask./mask); shading flat; 
    caxis([8 30]); 

    m_text(-96.5,43,extitle{ex},'fontsize',10,'fontweight','bold');
end

%=========================================================================%
%=                             Add Colorbar                              =%
%=========================================================================%
cax=axes('position',[0.7 0.05 0.1 0.65]); axis off
caxis([8.0 30])
hCbar     =colorbar('location','East','position',[0.7 0.05 0.1 0.67],...
                     'yaxislocation','right');
hCbarT=text(1.65,0.75,'dQ/dSST (Wm^{-2}K^{-1})','fontsize',12,...
                      'fontname','helvetica','fontweight','bold',...
                      'rotation',270);
set(hCbar,'linewidth',1,'xtick',[0:5:30],'ticklength',[.02 .02],...
          'fontsize',10,'fontname','helvetica','fontweight','bold');


%=========================================================================%
%=              Make set of text at top for model run info               =%
%=========================================================================%
tax=axes('position',[0.7 0.27 0.275 0.2],'xtick',[],'ytick',[]); box on;
set(tax,'position',[0.67 0.75 0.225 0.225]);
set(tax,'linewidth',1); 
hmodelname=text(0.03,0.9,'US East Coast ROMS','fontweight','bold',...
    'fontsize',12);
hNCSU     =text(0.065,0.75,'NC State Univeristy (2014)','fontsize',10);             
hmodelres =text(0.45,0.60,'\Delta x,\Delta y = 1/12^{o}','fontsize',10);          
hmodelvert=text(0.76,0.46,'k = 36','fontsize',10);
hmodeldate=text(0.975,0.2,'Oct-Dec dQ/dSST',...
                  'fontsize',14,'fontweight','bold',...
                  'horizontalalignment','right');

set(gcf,'color','w');
%axes(ax5); 
axes(ax2); axes(ax3); axes(ax1);
wysiwyg