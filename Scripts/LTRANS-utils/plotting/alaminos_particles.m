%============================ plot_site_locations.m ===========================%
%=                                                                            =%
%=  Written by Austin C. Todd, NCSU (2015) for personal use                   =%
%=                                                                            =%
%=  Program loads US East Coast ROMS Bathymetry and grid and plots the        =%
%=  position of particles at any given time. Frames are saved for making an   =%
%=  animation of particle movement for passive and bathymodiolus particles.   =%
%=                                                                            =%
%==============================================================================%

addpath('~/MYMATLAB/');
addpath(genpath('~/MYMATLAB/ROMS-matlab/'));
addpath(genpath('/Users/actodd/MYMATLAB/m_map'));
format long g; format compact;

%==============================================================================%
%=                     Set file input/output directories                      =%
%==============================================================================%
plotdir  ='/Volumes/Black_box/Data/PLOTS/Connectivity/frames/';
gridfile  = '/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc';
ltransdir='/Users/actodd/Documents/Work/Projects/Connectivity/Data/';
readdata=0;
site=1;

%==============================================================================%
%=                         Load ROMS grid information                         =%
%==============================================================================%
disp('Reading in ROMS output');

if readdata
ncid=netcdf.open(gridfile,'nowrite');
 lon        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'    ));
 lat        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'    ));
 mask       =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'   ));
 h          =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'          ));
netcdf.close(ncid);
end

%==============================================================================%
%=       Set the site lat/lon pairs (NOTE: sites 4-6 will not be used)        =%
%==============================================================================%
sitename='Alaminos Canyon';
sitelat=lat(47,210); % site1
sitelon=lon(47,210);
sitename='Florida Escarpment';
sitelat=lat(144,201); % site 3
sitelon=lon(144,201);

sitelat=[26.354722, 27.723056, 26.03,      11.233367, 10.327833,...
	 10.327833, 32.490833, 32.9741617, 37.0985,   38.166667];
sitelon=[-94.496667, -91.275,    -84.915,    -59.34595,  -58.88895, ...
         -58.88895,  -76.188333, -75.921667, -74.623167, -73.833333];
sitenames={'Alaminos Canyon';'Brine Pool';'Florida Escarpment';'El Pilar';...
           'Orenoque A';'Orenoque B ';'Blake Ridge';'Cape Fear';...
	   'Norfolk Canyon';'Baltimore Canyon'};

sitelat=sitelat([1:3 7:10]); sitelon=sitelon([1:3 7:10]); 
sitenames=sitenames([1:3 7:10]);

sitelat(1)=26.35473; sitelon(1)=-94.49690;
sitelat(3)=lat(144,201); sitelon(3)=lon(144,201);

sitelat=sitelat(site);
sitelon=sitelon(site);
ssite=sprintf('%1i',site);

%==============================================================================%
%=                           Read in LTRANS output                            =%
%==============================================================================%
if readdata
disp('Reading in LTRANS data');    
ncid=netcdf.open([ltransdir,'passive/site',ssite,'.nc'],'nowrite');
  passive.lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
  passive.lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
  passive.dob  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dob'));
  passive.age  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age'));
  passive.depth=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'depth'));
  passive.color=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'color'));
netcdf.close(ncid);
ncid=netcdf.open([ltransdir,'bathymodiolus/site',ssite,'.nc'],'nowrite');
  bathymod.lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
  bathymod.lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
  bathymod.dob  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dob'));
  bathymod.age  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age'));
  bathymod.depth=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'depth'));
  bathymod.color=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'color'));
netcdf.close(ncid);
end
%==============================================================================%
%=                       Plot grid and station locations                      =%
%==============================================================================%
for tdim=366:size(passive.lat,2)
  
  %============================================================================%
  %=                          Plot Passive Particles                          =%
  %============================================================================%
  clf;
  ax1=axes('position',[0.06 0.15 0.45 0.8]);
  m_proj('lambert','long',[-98 -80],...%max(max(lon(1:end-25,:)))],...
                 'lat' ,[18 31]); %max(max(lat(1:end-25,:)))]);
  m_pcolor(lon,lat,mask); shading flat; hold on;
  caxis([-3 1]); colormap(gray)
  m_grid('tickdir','out');
  
  %--- Bathymetry Contours ---%
  m_contour(lon,lat,h.*mask,[100 500 1000 2000 3000],'color',[.6 .6 .6],...
           'linewidth',1)
  m_contour(lon,lat,mask,[0 0],'k','linewidth',1)                                     

  %--- Make scatter of the particles ---%
  a=find(passive.dob<=(tdim-1)*86400); %<--- Particles released before tdim
  b=find(passive.age(a,tdim)<max(passive.age(:)));
  m_plot(passive.lon(b,tdim),passive.lat(b,tdim),'bo','markersize',4,...
         'markerfacecolor','b','markeredgecolor','b'); 
  
  %--- Site Locations ---%
  m_plot(sitelon,sitelat,'rs','markersize',5,'markerfacecolor','r',...
          'markeredgecolor','k')
  
  %--- Add title ---%
  title('Currents Only','fontsize',14,'fontweight','bold');
  
  %============================================================================%
  %=                      Plot Bathymodiolus Particles                        =%
  %============================================================================%  
  ax2=axes('position',[0.52 0.15 0.45 0.8]);
  m_proj('lambert','long',[-98 -80],...%max(max(lon(1:end-25,:)))],...
                 'lat' ,[18 31]); %max(max(lat(1:end-25,:)))]);
  m_pcolor(lon,lat,mask); shading flat; hold on;
  caxis([-3 1]); colormap(gray)
  m_grid('tickdir','out','yticklabel',[]);
  
  %--- Bathymetry Contours ---%
  m_contour(lon,lat,h.*mask,[100 500 1000 2000 3000],'color',[.6 .6 .6],...
           'linewidth',1)
  m_contour(lon,lat,mask,[0 0],'k','linewidth',1)                                     

  %--- Make scatter of the particles ---%
  a=find(bathymod.dob<=(tdim-1)*86400); %<--- Particles released before tdim
  b=find(bathymod.age(a,tdim)<max(bathymod.age(:)));
  m_plot(bathymod.lon(b,tdim),bathymod.lat(b,tdim),'mo','markersize',4,...
         'markerfacecolor','m','markeredgecolor','m'); 
  
  %--- Site Locations ---%
  m_plot(sitelon,sitelat,'rs','markersize',5,'markerfacecolor','r',...
          'markeredgecolor','k')
  
  %--- Add title ---%
  title('Bathymodiolus','fontsize',14,'fontweight','bold');

  %============================================================================%
  %=                      Mike Time Scale along bottom                        =%
  %============================================================================%  
  day_count=[0,31,28,31,30,31,30,31,31,30,31,30,31];
  months={'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};
  ax3=axes('position',[0.065 0.17 0.9 0.1]);
  
  %--- Major axes of line ---%
  plot([1 365],[0 0],'k','linewidth',2); hold on;
  plot([1 1],[-1 1],'k','linewidth',2) 
  plot([365 365],[-1 1],'k','linewidth',2)
   
  %--- Plot thick red line to indicate current date
  plot([1 1].*(mod(tdim-1,365)+1),[-1 1],'r','linewidth',8);

  %--- Ticks at Month transitions ---%
  for i=1:12
    plot([1 1].*sum(day_count(1:i+1)),[-.5 .5],'k','linewidth',1);
    text([sum(day_count(1:i))+day_count(i+1)/2],0.75,months{i},...
           'fontsize',12,'horizontalalignment','center');
  end
  text(mod(tdim-1,365)+1,-2,['Day ',sprintf('%03i',tdim)],'fontsize',12,...
       'horizontalalignment','center');%,'verticalalignment','center');
  axis off; hold on; axis([-1 367 -2 2]);
  
  
%  drawnow;
 % wysiwyg;
%  pause(0.25)
%return;
 print('-dpng','-r150','-painters',[plotdir,'frame_',sprintf('%03i',tdim),'.png']);

end