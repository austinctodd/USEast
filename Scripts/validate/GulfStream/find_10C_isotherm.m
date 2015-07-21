%========================= find_15C_isotherm.m ===========================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This reads the temperature field from the US East water age model    =%
%=  to find the location of the 10 degree isotherm at the surface (which =%
%=  is indicative of the Shelf Slope Front). The position is extracted   =%
%=  in the form of a latitude at each grid longitude (using interpolation=%
%=  routines.                                                            =%
%=                                                                       =%
%=========================================================================%

%--- Add various libraries and paths ---%
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/home/actodd/MYMATLAB/');
addpath(genpath('/home/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=========================================================================%
%=                  Set file input/output directories                    =%
%=========================================================================%
plotdir='/he_data/he/actodd/PLOTS/USeast/GS_position/';
datadir='/gpfs_common/he_share/actodd/USeast-age/output/clim/';
matdir ='/he_data/he/actodd/DATA/eke/';

%=========================================================================%
%=                      Load ROMS grid information                       =%
%=========================================================================%
disp('Reading in ROMS output');

ncid=netcdf.open([datadir,'year1/useast_his_0001.nc'],'nowrite');
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

%--- Set up depth lat/lon points ---%
lon2=lon(230:402,320:465);
lat2=lat(230:402,320:465);
h2  =h(  230:402,320:465);

%=========================================================================%
%=                       Loop through time record                        =%
%=========================================================================%

count=ones(1,173);

fname=[datadir,'daily_avg.nc'];
ncid=netcdf.open(fname,'nowrite');
for j=1:364

  disp(['Day ',sprintf('%03i',j)]);

  %--- Load temp and SSH, calculate depth ---%
  temp=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'temp'),...
                       [229 319 35 j-1],[173 146 1 1]);
  temp=squeeze(double(temp));


  %--- Go along longitude direction, find where 15C occurs ---% 
  for i=1:173
    a=find(isnan(temp)==0);
    if length(a)>0
      x=find( abs(temp(i,:)-10)==nanmin(abs(temp(i,:)-10)));
      templat(j,i)=lat2(i,x(1));
      tempval(j,i)=temp(i,x(1));
    end
  end
end
netcdf.close(ncid);

save /he_data/he/actodd/DATA/GS_position/slope_front_lats.mat templat tempval count -MAT

