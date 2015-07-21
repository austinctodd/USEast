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
datadir='/Volumes/Black_box/USeast-age/output/testdQdSST/';
climname='/Volumes/Black_box/Data/USeast/Data/clim/USeast-clim-monthly.nc';

expname={'ERA';'ERA-HYCOM';'CFSR';'CFSR-HYCOM'};

%=========================================================================%
%=                      Load ROMS grid information                       =%
%=========================================================================%
disp('Reading in ROMS output');

ncid=netcdf.open([datadir,expname{1},'/useast_his.nc'],'nowrite');
 Vstretching=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vstretching'));
 Vtransform =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vtransform' ));
 theta_s    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_s'    ));
 theta_b    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_b'    ));
 hc         =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hc'         ));
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
%=                           Load HYCOM temps                            =%
%=========================================================================%
ncid=netcdf.open(climname,'nowrite');
  hycom.time=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'));
  temp=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'temp'),...
                                          [0 0 35 0],[402-25 482 1 12]);
netcdf.close(ncid);
temp=squeeze(temp);

%--- Add in times for December before and January after ---%
temp(:,:,2:13)=temp(:,:,1:12); hycom.time(2:13)=hycom.time(1:12);
temp(:,:,1   )=temp(:,:,13);   hycom.time(1   )=hycom.time(13)-365.25;
temp(:,:,14  )=temp(:,:,2 );   hycom.time(14  )=hycom.time(2 )+365.25;
return;
%=========================================================================%
%=                   Interpolate HYCOM temps to daily                    =%
%=========================================================================%
for i=1:402-25
  for j=1:482
    hycom.temp(i,j,:)=interp1(hycom.time,squeeze(temp(i,j,:)),[1:366]);
  end
end
        
%=========================================================================%
%=                     Set different axes properties                     =%
%=========================================================================%
ax1pos(1,:)=[0.035 0.53 0.35 0.465];
ax1pos(2,:)=[0.340 0.53 0.35 0.465];
ax1pos(3,:)=[0.035 0.03 0.35 0.465];
ax1pos(4,:)=[0.340 0.03 0.35 0.465];
ax1pos(5,:)=[0.645 0.53 0.35 0.465];

%=========================================================================%
%=                          Load mean velocities                         =%
%=========================================================================%
for tdim=1:366
%  figure(1);
  clf; %colormap(cmap);
  mn=floor(tdim./31)+1;
  smn=sprintf('%1i',mn); 
  day=mod(tdim,31);
  disp(['Month ',smn,', Day ',sprintf('%02i',day)]);

  for ex=1:4

    %====================================================================%
    %=                   Load in surface velocities                     =%
    %====================================================================%
    
    ncid=netcdf.open([datadir,expname{ex},'/useast_his.nc'],'nowrite');
      temp=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'temp'),...
                                             [0 0 35 tdim-1],...
                                             [402-25 482 1 1]);
      temp=double(temp);
    netcdf.close(ncid);
  		   

    %====================================================================%
    %=                  Main plot panel for each exp                    =%
    %====================================================================%
    eval(['ax',sprintf('%1i',ex),'=axes(''position'',ax1pos(',...
               sprintf('%1i',ex),',:));']);
    m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                     'lat' ,[min(min(lat)) max(max(lat))]);
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
    pcc=m_pcolor(lon,lat,temp.*mask./mask); shading flat; 
    caxis([0.0 30]); %alpha(pcc,0.65);

    m_text(-96.5,43,expname{ex},'fontsize',10,'fontweight','bold');
%    m_text(-96.5,38,['Tnu2=',sprintf('%02i',tnu2(1))],'fontsize',10,...
%                                                     'fontweight','bold');
%    m_text(-96.5,36,['Vis2=',sprintf('%02i',vis2  )],'fontsize',10,...
%                                                     'fontweight','bold');

  end

  %======================================================================%
  %=                    Main plot panel for HYCOM                       =%
  %======================================================================%
  ax5=axes('position',ax1pos(5,:));
  m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                   'lat' ,[min(min(lat)) max(max(lat))]);
  m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
         'xticklabels',char({'96^{o}W';'88^{o}W';'80^{o}W';'72^{o}W';'64^{o}W'}),...
         'yticklabels',char({'';'';'32^{o}N ';'40^{o}N '}),'fontsize',10)
  hold on;
  m_contour(lon,lat,mask,[0 0],'k','linewidth',1.5);
  pcc=m_pcolor(lon,lat,squeeze(hycom.temp(i,j,tdim)).*mask./mask);
  shading flat; 
  caxis([0.0 30]); 
  m_text(-96.5,43,['HYCOM'],'fontsize',10,'fontweight','bold');


  %--- Add Colorbar (scaled by current day) ---%
  cax=axes('position',[0.7 0.08 0.25 0.35]); axis off
  %[0.82 0.075 0.075 0.7]); axis off;
  caxis([0.0 30])
  hCbar     =colorbar('location','South','position',[0.7 0.175 0.275 0.08],...
                      'XAxisLocation','Bottom');
  hCbarT=text(0.135,0.1,'Sfc Temperature (C)','fontsize',12,...
                                       'fontname','helvetica',...
                                       'fontweight','bold');
  set(hCbar,'linewidth',1,'xtick',[0:5:30],'ticklength',[.02 .02],...
            'fontsize',10,'fontname','helvetica','fontweight','bold');


  %--- Make set of text at top for model run info ---%
  tax=axes('position',[0.7 0.27 0.275 0.2],'xtick',[],'ytick',[]); box on;
  set(tax,'linewidth',1); 
  hmodelname=text(0.2,0.9,'US East Coast ROMS',...
                              'fontweight','bold','fontsize',12);
  hNCSU     =text(0.225,0.75,'NC State Univeristy (2014)','fontsize',10);             
  hmodelres =text(0.56,0.60,'\Delta x,\Delta y = 1/12^{o}','fontsize',10);          
  hmodelvert=text(0.79,0.46,'k = 36','fontsize',10);
  hmodeldate=text(0.975,0.2,[datestr(datenum(2013,1,tdim))],...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');

  set(gcf,'color','w');
  axes(ax5); axes(ax2); axes(ax3); axes(ax1);
  wysiwyg
  eval(['export_fig ',plotdir,'temp.',sprintf('%03i',tdim),'.png ',...
       '-png -r150 -painters']);
end

