%============================= ROMS_EKE.m ================================%
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
datadir='/he_data/he/actodd/DATA/eke/'; % For saving eke.mat files

disp('Reading in ROMS output');
pdir='/he_data/he/zhigang/Project/USeast/Data/';

ncid=netcdf.open([pdir,'USeast-clm.nc'],'nowrite');
lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho' ));
%lonu =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_u'   ));
%lonv =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_v'   ));
lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho' ));
%latu =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_u'   ));
%latv =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_v'   ));
%mask =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'));
%masku=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_u'  ));
%maskv=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_v'  ));
h    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'       ));

%=========================================================================%
%=                       Find mean u and v fields                        =%
%=========================================================================%
disp('Calculating mean velocity fields');
meanu=zeros(36,401,482);  meanv=zeros(36,402,481);

for j=186-30:186
  %--- Read u and add to mean ---%
  u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 0 j-1],[401 482 36 1]);
  meanu(:,:,:)=meanu(:,:,:)+double(permute(u,[3 1 2]));

  %--- Read v and add to mean ---%
  u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 0 j-1],[402 481 36 1]);
  meanv(:,:,:)=meanv(:,:,:)+double(permute(u,[3 1 2]));

end
%meanu=meanu./186; meanv=meanv./186;
meanu=meanu./31; meanv=meanv./31;

eval(['save ',datadir,'hycom_JUNE_meanvls.mat meanu meanv -MAT']);
%load([datadir,'mean_vels.mat']);

%=========================================================================%
%=                       Loop through to calculate EKE                   =%
%=========================================================================%
disp('Calculating EKE');

ke=zeros(36,402,482); count=0;
for i=6:6
  for j=1:31

    %--- Read u and subtract mean ---%
    u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 0 j-1],[401 482 36 1]);
    u=double(permute(u,[3 1 2]))-meanu(:,:,:);
     
    %--- Read v and subtract mean ---%
    v=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 0 j-1],[402 481 36 1]);
    v=double(permute(v,[3 1 2]))-meanv(:,:,:);
      
    for k=1:36
      %--- Regrid velocities to rho-grid ---%
      [uu,vv]=regridromsvels2d(squeeze(u(k,:,:))',squeeze(v(k,:,:))');
      uu=uu'.*100; vv=vv'.*100;
      temp(1,:,:)=(0.5*( uu.^2 + vv.^2 ));
  
      %--- Calculate EKE ---%
      ke(k,:,:)=ke(k,:,:)+temp; clear temp
    end
    count=count+1;      
  end
end
ke=ke./count;
eval(['save ',datadir,'eke_JUNE_hycom.mat ke -MAT']);


