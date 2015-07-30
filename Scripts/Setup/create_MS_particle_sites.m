%=========================================================================%
%
%  Program:  create_ltrans_MS_release.m
%  Author: Austin C Todd (NCSU) 2015
%
%  Method:  This program creates particle release times and locations for
%           six different sites near the Mississippi and Atchafalaya River
%           mouths in the US East model for particle releases.
%
%=========================================================================%

outputdir = '/Volumes/Black_box/Data/LTRANS/input/USEast/Mississippi/';
grid_file = '/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc';

%------------------------------------------------
% Open particle release files
%------------------------------------------------
fid1=fopen([outputdir,'AR.csv'],'w');
fid2=fopen([outputdir,'MR1.csv'],'w');
fid3=fopen([outputdir,'MR2.csv'],'w');
fid4=fopen([outputdir,'MR3.csv'],'w');
fid5=fopen([outputdir,'MR4.csv'],'w');
fid6=fopen([outputdir,'MR5.csv'],'w');

%------------------------------------------------
% Open ROMS Grid file
%------------------------------------------------
ncid=netcdf.open(grid_file,'nowrite');
  lon  = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'));
  lat  = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'));
  lonu = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_u'));
  latu = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_u'));
  lonv = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_v'));
  latv = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_v'));
  h    = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'));
  mask = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'));
netcdf.close(ncid);

%------------------------------------------------
% Define particle release locationss
%------------------------------------------------
plon(1)=lon( 79,243); plat(1)=lat(79,243);  
plon(2)=lon( 99,237); plat(2)=lat(99,237);
plon(3)=lon(101,237); plat(3)=lat(101,237);
plon(4)=lon(102,242); plat(4)=lat(102,242);
plon(5)=lon(103,238); plat(5)=lat(103,238);
plon(6)=lon(104,239); plat(6)=lat(104,239);

%------------------------------------------------
% Output to file
%------------------------------------------------
for i=0:1/8:365.25
  for j=1:10
    fprintf(fid1,'%8.5f %9.5f -0.50 %d\n',[plat(1) plon(1) i*86400]);
    fprintf(fid2,'%8.5f %9.5f -0.50 %d\n',[plat(2) plon(2) i*86400]);
    fprintf(fid3,'%8.5f %9.5f -0.50 %d\n',[plat(3) plon(3) i*86400]);
    fprintf(fid4,'%8.5f %9.5f -0.50 %d\n',[plat(4) plon(4) i*86400]);
    fprintf(fid5,'%8.5f %9.5f -0.50 %d\n',[plat(5) plon(5) i*86400]);
    fprintf(fid6,'%8.5f %9.5f -0.50 %d\n',[plat(6) plon(6) i*86400]);
  end
end

%------------------------------------------------
% Close files
%------------------------------------------------
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);
