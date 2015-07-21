%============================== AVISO_EKE.m ==============================%
%=                                                                       =%
%=  Written by Austin C. Todd, NCSU (2015) for personal use              =%
%=                                                                       =%
%=  Program loads AVISO Geostrophic velocity anomolies and calculates    =%
%=  the EKE over the US East ROMS model domain.                          =%
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
           'dataset-duacs-dt-global-allsat-msla-uv?u'];
latinds  ='[383:1:548]';
loninds  ='[1043:1:1204]';

%=========================================================================%
%=                      Load ROMS lat/lon from file                      =%
%=========================================================================%
fname='/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc';
ncid=netcdf.open(fname,'nowrite');
  roms.lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho' ));
  roms.lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho' ));
  roms.mask=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'));
netcdf.close(ncid)

%=========================================================================%
%=                           Loop through files                          =%
%=========================================================================%
EKE=zeros(166,162);
count=0;

for tdim = 0:7669

  stdim=sprintf('%i',tdim);
  disp(['Time ',sprintf('%04i',tdim+1),'/7670']);
  
  %--- Construct url and load data ---$  
  fname=[AVISOurl,'[',stdim,':1:',stdim,']',latinds,loninds,',v',...
                  '[',stdim,':1:',stdim,']',latinds,loninds,',crs'];
  loaddap(fname);
     
  %--- Data management ---%
  u.u(find(u.u==crs))=nan;
  u.u=u.u./0.0001./crs;
  
  v.v(find(v.v==crs))=nan;
  v.v=v.v./0.0001./crs;

  %--- Calculate EKE ---% 
  EKE = EKE + 0.5*(u.u.^2 + v.v.^2);
  count=count+1;
end


