function plot_procedure(lon,lat,mask,var,C,textlab,textlab2)
%=========================================================================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This routine is the core of the plotting procedure for the KE and    =%
%=  EKE plotting.  It uses m_map.                                        =%
%=                                                                       =%
%=========================================================================%

addpath('/home/actodd/MYMATLAB/');
addpath(genpath('/home/actodd/MYMATLAB/m_map'));

%--- Plot the KE ---%
figure(1); clf;
ax1=axes('position',[0.01 0.05 0.8 0.93]);
m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                 'lat' ,[min(min(lat)) max(max(lat))]);
m_grid('tickdir','on');
hold on;
m_pcolor(lon,lat,var.*mask./mask); shading flat;
m_contour(lon,lat,mask,[0 0],'k');
m_text(-95.5,38,textlab2,'fontsize',14,'fontweight','bold');
m_text(-92.5,36.75,'(surface)','fontsize',14,'fontweight','bold');
caxis(C);

%--- Add Colorbar (scaled by current day) ---%
cax=axes('position',[0.82 0.075 0.075 0.7]); axis off;
caxis(C);
hCbar     =colorbar('position',[0.82 0.075 0.075 0.7]);
hCbarT=text(-0.5,-0.05,'log(KE) (cm^{2}s^{-2})','fontsize',12,...
                                     'fontname','helvetica',...
                                     'fontweight','bold');
set(hCbar,'linewidth',1,'ytick',[min(C):(max(C)-min(C))./10:max(C)],...
          'fontsize',10,'fontname','helvetica','fontweight','bold');


%--- Make set of text at top for model run info ---%
tax=axes('position',[0.7 0.8 0.275 0.175]); axis off;
hmodelname=text(0.165,0.95,'US East Coast ROMS',...
                            'fontweight','bold','fontsize',12);
hNCSU     =text(0.18,0.775,'NC State Univeristy (2014)','fontsize',10);             
hmodelres =text(0.51,0.56,'\Delta x,\Delta y = 1/10^{o}','fontsize',10);          
hmodelvert=text(0.75,0.395,'k = 36','fontsize',10);
hmodeldate=text(0.925,0.2,textlab,...
                    'fontsize',14,'fontweight','bold',...
                    'horizontalalignment','right');

set(gcf,'color','w');
wysiwyg;
return;
