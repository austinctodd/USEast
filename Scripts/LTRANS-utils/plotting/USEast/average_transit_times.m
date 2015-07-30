%=========================================================================%
%
% PROGRAM: average_transit_times.m
%
% PURPOSE: Program reads in LTRANS output and calculates the length of time
%          that each particle spends within the specified residence time
%          domain. 
%
% AUTHOR: Austin C. Todd (NCSU), 9 March 2015
%
%=========================================================================%

%=========================================================================%
% Set file directories
%=========================================================================%
LTRANS_dir ='/Volumes/Black_box/Data/LTRANS/output/Mississippi/';
ROMS_dir   ='/Volumes/Black_box/Data/USeast/Data/grd/';
PLOT_dir   ='/Volumes/Black_box/Data/PLOTS/LTRANS/Mississippi/';

%=========================================================================%
% Create file names for particle files and grid file
%=========================================================================%
AR.fname  =[LTRANS_dir,'AR_v2.txt'];
MS1.fname =[LTRANS_dir,'MS1_v2.txt'];
MS2.fname =[LTRANS_dir,'MS2_v2.txt'];
MS3.fname =[LTRANS_dir,'MS3_v2.txt'];
MS4.fname =[LTRANS_dir,'MS4_v2.txt'];
MS5.fname =[LTRANS_dir,'MS5_v2.txt'];
grid_file =[ROMS_dir,'grid_GOM_shelf_scope.nc'];

% Rivers
rivers={'AR';'MS1';'MS2';'MS3';'MS4';'MS5'}

%=========================================================================%
% Open grid file and read lat/lon/scope
%=========================================================================%
disp(['Reading data from ROMS file']);
ncid=netcdf.open(grid_file,'nowrite');
  lon   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'));
  lat   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'));
  scope = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'scope_rho'));
  mask  = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'));
netcdf.close(ncid);

%=========================================================================%
% Loop through each river mouth to read in and average data
%=========================================================================%
trans_time=zeros(200,161);
trans_part=zeros(200,161);
expos_time=zeros(200,161);
expos_part=zeros(200,161);

for f=1:6
  if f==1
    fname=AR.fname;
  else
    eval(['fname=[',rivers{f},'.fname];']);
  end
  disp(['Reading data from ',fname]);
  
  [ii,jj,tt,tp,et,ep]=textread(fname,'%n %n %f %n %f %n\n');
  for i=1:32200
    trans_time(ii(i),jj(i)-99)=trans_time(ii(i),jj(i)-99)+tt(i);
    trans_part(ii(i),jj(i)-99)=trans_part(ii(i),jj(i)-99)+tp(i);
    expos_time(ii(i),jj(i)-99)=expos_time(ii(i),jj(i)-99)+et(i);
    expos_part(ii(i),jj(i)-99)=expos_part(ii(i),jj(i)-99)+ep(i);
  end
end

