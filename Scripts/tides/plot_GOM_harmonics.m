%========================= plot_GOM_harmonics.m ==========================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  Program plots the tidal phases and amplitudes on the same plot with  =%
%=  similar colors and domain size as seen in figure 5 of Guoillon et al =%
%=  (2010). Requires the tidal harmonics to have been computed using     =%
%=  t_tide routines.                                                     =%
%=                                                                       =%
%=========================================================================%

%=========================================================================%
%=                        Load dependent libraries                       =%
%=========================================================================%
addpath(genpath('/he_data/he/zhigang/Matlab_Codes/mexcdf2009/'));
addpath('/home/actodd/MYMATLAB/');

%=========================================================================%
%=                      Set filenames and load data                      =%
%=========================================================================%
exp='LTP';  coding='newcode_orig';

dir=['/gpfs_share/rhe/actodd/USeast-tide/output/',exp,'/',coding,'/'];
load([dir,'GOM_harmonic_analysis.mat']);
load /home/actodd/ROMS-utils/tides/flav_cmap.mat

fname=[dir,'useast_his.nc'];
if (size(O1phase,1)<400)
  lon =nc_varget(fname,'lon_rho' ,[80 9],[183 177]); 
  lat =nc_varget(fname,'lat_rho' ,[80 9],[183 177]);
  mask=nc_varget(fname,'mask_rho',[80 9],[183 177]);
else
  lon =nc_varget(fname,'lon_rho' ); 
  lat =nc_varget(fname,'lat_rho' );
  mask=nc_varget(fname,'mask_rho');
end

%=========================================================================%
%=                     Create plots and save to file                     =%
%=========================================================================%

%--- O1 phases and amplitudes in the GOM ---%
figure(1); clf
pcolor(lon,lat,O1phase); shading flat; axis xy; axis image
axis([-98 -81 15.5 31])
caxis([0 360])
colormap(cmap)
hold on
[co1,ho1]=contour(lon,lat,O1amp,[0.05:0.05:0.3],'k');
set(gca,'tickdir','out','fontsize',10,'linewidth',1)
clabel(co1,ho1,'fontsize',10)
colorbar
title('O1 phases and amplitudes','fontsize',12,'fontweight','bold')
wysiwyg
set(gcf,'Color','w');
eval(['export_fig /home/actodd/PLOTS/',exp,'/',coding,'/O1_GOM.png ',...
      '-png -r150 -painters']);

%--- K1 phases and amplitudes in the GOM ---%
figure(2); clf
pcolor(lon,lat,K1phase); shading flat; axis xy; axis image
axis([-98 -81 15.5 31])
caxis([0 360])
colormap(cmap)
hold on
[ck1,hk1]=contour(lon,lat,K1amp,[0.05:0.05:0.3],'k');
set(gca,'tickdir','out','fontsize',10,'linewidth',1)
clabel(ck1,hk1,'fontsize',10)
colorbar
title('K1 phases and amplitudes','fontsize',12,'fontweight','bold')
wysiwyg
set(gcf,'Color','w');
%export_fig /home/actodd/PLOTS/LTP/K1_GOM.png -png -r150 -painters
eval(['export_fig /home/actodd/PLOTS/',exp,'/',coding,'/K1_GOM.png ',...
      '-png -r150 -painters']);


%--- M2 phases and amplitudes in the GOM ---%
figure(3); clf
pcolor(lon,lat,M2phase); shading flat; axis xy; axis image
axis([-98 -81 15.5 31])
caxis([0 360])
colormap(cmap)
hold on
[cm2,hm2]=contour(lon,lat,M2amp,[0.05:0.05:0.3],'k');
set(gca,'tickdir','out','fontsize',10,'linewidth',1)
clabel(cm2,hm2,'fontsize',10)
colorbar
title('M2 phases and amplitudes','fontsize',12,'fontweight','bold')
wysiwyg
set(gcf,'Color','w');
%export_fig /home/actodd/PLOTS/LTP/M2_GOM.png -png -r150 -painters
eval(['export_fig /home/actodd/PLOTS/',exp,'/',coding,'/M2_GOM.png ',...
      '-png -r150 -painters']);


%--- S2 phases and amplitudes in the GOM ---%
figure(4); clf
pcolor(lon,lat,S2phase); shading flat; axis xy; axis image
axis([-98 -81 15.5 31])
caxis([0 360])
colormap(cmap)
hold on
[cs2,hs2]=contour(lon,lat,S2amp,[0.05:0.05:0.3],'k');
set(gca,'tickdir','out','fontsize',10,'linewidth',1)
clabel(cs2,hs2,'fontsize',10)
colorbar
title('S2 phases and amplitudes','fontsize',12,'fontweight','bold')
wysiwyg
set(gcf,'Color','w');
%export_fig /home/actodd/PLOTS/LTP/S2_GOM.png -png -r150 -painters
eval(['export_fig /home/actodd/PLOTS/',exp,'/',coding,'/S2_GOM.png ',...
      '-png -r150 -painters']);

