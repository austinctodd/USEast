%=========================== find_LC_contour.m ===========================%
%=                                                                       =%
%=  Written by Austin C. Todd, NCSU (2015) for personal use              =%
%=                                                                       =%
%=  Program loads US East Coast ROMS Bathymetry and grid and plots the   =%
%=  locations that will be used for particle advection experiments.      =%
%=                                                                       =%
%=========================================================================%

addpath(genpath('/Users/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/Users/actodd/MYMATLAB/');
addpath(genpath('/Users/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=========================================================================%
%=                   Set file input/output directories                   =%
%=========================================================================%
AVISOdir ='/Volumes/Black_box/Data/AVISO/SLA/';
ROMSdir  ='/Volumes/Black_box/Data/USeast-age/output/clim/averages/';

%=========================================================================%
%=                   Set file input/output directories                   =%
%=========================================================================%
AVISOurl =['http://opendap.aviso.altimetry.fr/thredds/dodsC/',...
           'dataset-duacs-dt-global-allsat-msla-h?sla'];
latinds  ='[430:1:485]';
loninds  ='[1045:1:1120]';

%=========================================================================%
%=                       Load ROMS mean zeta field                       =%
%=========================================================================%
ncid=netcdf.open([ROMSdir,'mean_vels.nc'],'nowrite');
  roms.zeta=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mean_zeta'),[0 104],[198 163]);
netcdf.close(ncid)
roms.zeta=double(roms.zeta');

%=========================================================================%
%=                      Load ROMS lat/lon from file                      =%
%=========================================================================%
fname='/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc';
ncid=netcdf.open(fname,'nowrite');
  roms.lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho' ),[0 104],[1 163]);
  roms.lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho' ),[0   0],[198 1]);
  roms.mask=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'),[0 104],[198 163]);
netcdf.close(ncid)

%=========================================================================%
%=                     Load AVISO lat/lons from file                     =%
%=========================================================================%
fname=[AVISOdir,'1993/dt_global_allsat_msla_h_19931210_20140106.nc'];
ncid=netcdf.open(fname,'nowrite');
  aviso.lat=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'),430,56);
  aviso.lon=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'),1045,76);
netcdf.close(ncid)
[x,y]=meshgrid(double(aviso.lon-360),double(aviso.lat));

%=========================================================================%
%=                           Loop through files                          =%
%=========================================================================%
LCcnt=zeros(length(roms.lon),length(roms.lat));
count=0;

for tdim = 0:7669

  stdim=sprintf('%i',tdim);
  disp(['Time ',sprintf('%04i',tdim+1),'/7670']);
  
  %--- Construct url and load data ---$  
  fname=[AVISOurl,'[',stdim,':1:',stdim,']',latinds,loninds];
  loaddap(fname);
  aviso.sla=sla.sla;
  aviso.time(tdim+1)=sla.time;
  clear sla;
      
  %--- Data management ---%
  aviso.sla=double(aviso.sla).*.0001;
  aviso.sla(find(aviso.sla<=-21474))=nan;

  %--- Interpolate AVISO SLA to ROMS Grid ---%
  sla=griddata(x,y,aviso.sla,roms.lon,roms.lat);
  sla=(sla+roms.zeta).*roms.mask'./roms.mask';
      
  %--- Contour to obtain Loop Current & Ring positions --%
  clf; [c,h]=contour(roms.lon,roms.lat,sla,[0.17 0.17]);
      
  %--- Loop through contour points and add 
  LCpts=zeros(length(roms.lon),length(roms.lat));
  for i=1:length(c)
    if (c(2,i) > 18.0 && c(2,i)<roms.lat(end) &&...
        c(1,i) <-79.8 && c(1,i)>roms.lon(1  ))
           platind=max(find(roms.lat<=c(2,i)));
           plonind=floor((c(1,i)-roms.lon(1))*10.2790697675)+1;

       LCpts(plonind,platind)=1;
    end
  end
  LCcnt=LCcnt+LCpts;
  count=count+1;
end


