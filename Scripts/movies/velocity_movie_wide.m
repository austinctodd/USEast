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
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/home/actodd/MYMATLAB/');
addpath(genpath('/home/actodd/MYMATLAB/m_map'));
format long g; format compact;
set(0,'defaultfigurerenderer','zbuffer');

%=========================================================================%
%=                  Set file input/output directories                    =%
%=========================================================================%
plotdir='/he_data/he/actodd/PLOTS/USeast/velocities/sfc/frames/';
datadir='/he_data/he/zhigang/Project/USeast/Out_exp';
matdir ='/he_data/he/actodd/DATA/eke/';

%=========================================================================%
%=              Initialize constants for draw_arrow routine              =%
%=========================================================================%
sclv= 2;           % scaling= # degs / (1 cm/s) of velocity
col= [0.1 0.1 0.1];   % arrow color
lwd= 0.5;            % line width
beta= 20;             % angle between head and arrow
cf=0.4;

%========================================================================%
%=             Load in colormap used for velocity field                 =%
%========================================================================%
cmap=load('~/MYMATLAB/colormaps/density_colormap.txt');
cmap(end,:)=cmap(end-1,:);
for i=29:-1:1
  cmap(i,1)=cmap(i+1,1)+(1-cmap(30,1))/30;
  cmap(i,2)=cmap(i+1,2)+(1/30);
  cmap(i,3)=cmap(i+1,3)+(1-cmap(30,3))/30;
end

%=========================================================================%
%=                      Load ROMS grid information                       =%
%=========================================================================%
disp('Reading in ROMS output');

ncid=netcdf.open([datadir,'22/his_0001.nc'],'nowrite');
 Vstretching=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vstretching'));
 Vtransform =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vtransform' ));
 theta_s    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_s'    ));
 theta_b    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_b'    ));
 hc         =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hc'         ));
 lon        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'    ));
 lat        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'    ));
 mask       =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'   ));
 h          =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'          ));
netcdf.close(ncid);

[sizelt,sizeln]=size(mask);
nanmask=mask./mask;
zw=set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,36,5,h,h.*0,0);
h_layer=zw(:,:,2:end)-zw(:,:,1:end-1);

h100mask=mask; h50mask=mask;
h100mask(find(h<100))=0;
h50mask( find(h<50 ))=0;

%=========================================================================%
%=                     Set different axes properties                     =%
%=========================================================================%
ax1pos(1,:)=[0.035 0.53 0.35 0.465];
ax1pos(2,:)=[0.340 0.53 0.35 0.465];
ax1pos(3,:)=[0.035 0.08 0.35 0.465];
ax1pos(4,:)=[0.340 0.08 0.35 0.465];
ax1pos(5,:)=[0.645 0.53 0.35 0.465];

%=========================================================================%
%=                          Load mean velocities                         =%
%=========================================================================%
projs=[1,2,3,4,5,6,11,12,14,16,17,18,19,20,21,22];
pp=[7,8,9,16];

for tdim=1:31*6
  figure(1); clf; colormap(cmap);
  set(gcf,'position',[10 10 800 600])

  mn=floor(tdim./31)+1;
  smn=sprintf('%1i',mn); 
  day=mod(tdim,31)+1;
  disp(['Month ',smn,', Day ',sprintf('%02i',day)]);

  for ex=1:4

    %====================================================================%
    %=                   Load in surface velocities                     =%
    %====================================================================%
    fname=[datadir,sprintf('%02i',projs(pp(ex))),'/his_000',smn,'.nc'];
		   	   
    ncid=netcdf.open(fname,'nowrite');
      u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 35 mod(tdim,31)],...
                                                     [401 482 1 1]);
      v=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 35 mod(tdim,31)],...
                                                     [402 481 1 1]);
      tnu2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'nl_tnu2'));
      vis2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'nl_visc2'));
    netcdf.close(ncid);
    u(find(u>1e30))=nan; v(find(v>1e30))=nan;
  		   
    [uu,vv]=regridromsvels2d(double(u'),double(v'));
    uu=uu'; vv=vv';

    %====================================================================%
    %=                  Main plot panel for each exp                    =%
    %====================================================================%
    eval(['ax',sprintf('%1i',ex),'=axes(''position'',ax1pos(',...
               sprintf('%1i',ex),',:));']);
    m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                     'lat' ,[min(min(lat)) max(max(lat))]);
%    m_grid('tickdir','on');
    if (ex==1 || ex==3)
      m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
             'xticklabels',char({'96^{o}W';'88^{o}W';'80^{o}W';'72^{o}W';'64^{o}W'}),...
	     'yticklabels',char({'16^{o}N ';'24^{o}N ';'32^{o}N ';'40^{o}N ';}))
    else
      m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
             'xticklabels',char({'';'88^{o}W';'80^{o}W';'72^{o}W';'64^{o}W'}),...
	     'yticklabels','')
    end
    hold on;
    m_contour(lon,lat,mask,[0 0],'k','linewidth',1.5);
    pcc=m_pcolor(lon,lat,sqrt(uu.^2+vv.^2).*mask./mask); shading flat; 
    caxis([0.02 1.5]); alpha(pcc,0.65);

 
    % Add Scale Arrow
%    m_draw_arrow(-92,-92+1*sclv,41.75,41.75,cf,beta,col,lwd);
%    m_text(-96.5,43,'Surface Velocities','fontsize',10,'fontweight','bold');
%    m_text(-90  ,41.75,' 1 m/s','fontsize',10,'fontweight','bold');  
    m_text(-96.5,38,['Tnu2=',sprintf('%02i',tnu2(1))],'fontsize',10,...
                                                     'fontweight','bold');
    m_text(-96.5,36,['Vis2=',sprintf('%02i',vis2  )],'fontsize',10,...
                                                     'fontweight','bold');

  end

  %=======================================================================%
  %=                         Plot HYCOM velocities                       =%
  %=======================================================================%
  fname=[datadir(1:end-7),'/Data/USeast-clm.nc'];
		   	   
  ncid=netcdf.open(fname,'nowrite');
    u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 35 tdim-1],...
                                                   [401 482 1 1]);
    v=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 35 tdim-1],...
                                                   [402 481 1 1]);
  netcdf.close(ncid);
  u(find(u>1e30))=nan; v(find(v>1e30))=nan;
  	   
  [uu,vv]=regridromsvels2d(double(u'),double(v'));
  uu=uu'; vv=vv';

  %======================================================================%
  %=                    Main plot panel for HYCOM                       =%
  %======================================================================%
  ax5=axes('position',ax1pos(5,:));
  m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                   'lat' ,[min(min(lat)) max(max(lat))]);
  m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
         'xticklabels',char({'';'88^{o}W';'80^{o}W';'72^{o}W';'64^{o}W'}),...
         'yticklabels','')
  hold on;
  m_contour(lon,lat,mask,[0 0],'k','linewidth',1.5);
  pcc=m_pcolor(lon,lat,sqrt(uu.^2+vv.^2).*mask./mask); shading flat; 
  caxis([0.02 1.5]); alpha(pcc,0.65);
  %m_text(-96.5,43,'Surface Velocities','fontsize',10,'fontweight','bold');
  m_text(-96.5,38,['HYCOM'],'fontsize',10,'fontweight','bold');
  
  %--- Add Colorbar (scaled by current day) ---%
  cax=axes('position',[0.7 0.12 0.25 0.35]); axis off
  %[0.82 0.075 0.075 0.7]); axis off;
  caxis([0.02 1.5])
  hCbar     =colorbar('location','South','position',[0.7 0.215 0.275 0.08],...
                      'XAxisLocation','Bottom');
  hCbarT=text(0.135,0.1,'Current Speed (m/s)','fontsize',12,...
                                       'fontname','helvetica',...
                                       'fontweight','bold');
  set(hCbar,'linewidth',1,'xtick',[0.3:0.3:1.5],'ticklength',[.02 .02],...
            'fontsize',10,'fontname','helvetica','fontweight','bold');


  %--- Make set of text at top for model run info ---%
  tax=axes('position',[0.7 0.31 0.275 0.2],'xtick',[],'ytick',[]); box on;
  set(tax,'linewidth',1); 
  hmodelname=text(0.14,0.9,'US East Coast ROMS',...
                              'fontweight','bold','fontsize',12);
  hNCSU     =text(0.16,0.75,'NC State Univeristy (2014)','fontsize',10);             
  hmodelres =text(0.53,0.60,'\Delta x,\Delta y = 1/12^{o}','fontsize',10);          
  hmodelvert=text(0.78,0.46,'k = 36','fontsize',10);
  hmodeldate=text(0.975,0.2,[datestr(datenum(2013,1,tdim))],...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');

  set(gcf,'color','w');
  axes(ax5); axes(ax2); axes(ax3); axes(ax1);
  wysiwyg;

  eval(['export_fig ',plotdir,'velframe',sprintf('%02i',tdim),'.png ',...
       '-png -r150 -painters']);
end

