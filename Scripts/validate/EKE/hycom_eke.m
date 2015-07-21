%============================= HYCOM_EKE.m ===============================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program is designed to read in results from a US East Coast     =%
%=  ROMS simulation and calculate the KE and EKE. For the KE calculation =%
%=  the velocities are broken down into u=u_ + u', where u_ is the time  =%
%=  mean velocity, and u' is the time-varying perturbation velocity. The =%
%=  KE is then KE=0.5*(u_^2 + v_^2) and EKE=0.5*(u'^2 + v'^2).           =%
%=                                                                       =%
%=========================================================================%

%--- Add various libraries and paths ---%
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/home/actodd/MYMATLAB/');

%=========================================================================%
%=                      Load ROMS grid information                       =%
%=========================================================================%
disp('Reading in ROMS output');
pdir='/he_data/he/zhigang/Project/USeast/Data/';

ncid=netcdf.open([pdir,'USeast-clm.nc'],'nowrite');
lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho' ));
lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho' ));

%=========================================================================%
%=                       Find mean u and v fields                        =%
%=========================================================================%
disp('Calculating mean velocity fields');
meanu=zeros(401,482);  meanv=zeros(402,481);

for j=1:186
  %--- Read u and add to mean ---%
  u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 35 j-1],[401 482 1 1]);
  meanu(:,:)=meanu(:,:)+double(u);
  
  %--- Read v and add to mean ---%
  u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 35 j-1],[402 481 1 1]);
  meanv(:,:)=meanv(:,:)+double(u);
  
end
meanu=meanu./186; meanv=meanv./186;

save hycom_meanvls.mat meanu meanv -MAT

%=========================================================================%
%=                       Loop through to calculate EKE                   =%
%=========================================================================%
disp('Calculating EKE');

ke=zeros(186,402,482); count=1;

for j=1:186

  %--- Read u and subtract mean ---%
  u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 35 j-1],[401 482 1 1]);
  u=double(u)-meanu(:,:);
     
  %--- Read v and subtract mean ---%
  v=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 35 j-1],[402 481 1 1]);
  v=double(v)-meanv(:,:);
      
  %--- Regrid velocities to rho-grid ---%
  [uu,vv]=regridromsvels2d(u',v'); uu=uu'.*100; vv=vv'.*100;

  %--- Calculate EKE ---%
  ke(count,:,:)=0.5*( uu.^2 + vv.^2 );
  count=count+1;      

end
netcdf.close(ncid);

save eke_hycom.mat ke -MAT





