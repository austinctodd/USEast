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
pdir='/he_data/he/zhigang/Project/USeast/Out_exp';

ncid=netcdf.open([pdir,'01/his_0001.nc'],'nowrite');
lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho' ));
lonu =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_u'   ));
lonv =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_v'   ));
lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho' ));
latu =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_u'   ));
latv =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_v'   ));
mask =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'));
masku=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_u'  ));
maskv=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_v'  ));
h    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'       ));
netcdf.close(ncid);

%=========================================================================%
%=                       Find mean u and v fields                        =%
%=========================================================================%
disp('Calculating mean velocity fields');
%meanu=zeros(4,401,482);  meanv=zeros(4,402,481);

%load mean_vels.mat
%uu=meanu; meanu=zeros(6,401,482); meanu(1:5,:,:)=uu;
%uu=meanv; meanv=zeros(6,402,481); meanv(1:5,:,:)=uu;
%clear uu;

projs=[1,2,3,4,5,6,11,12,14,16,17,18,19,20,21,22];
meanu=zeros(length(projs),36,401,482);
meanv=zeros(length(projs),36,402,481);

eval(['load ',datadir,'mean_JUNE_vels.mat'])
meanu(length(projs),:,:,:)=0;
meanv(length(projs),:,:,:)=0;


for ex=length(projs):length(projs)
  disp(['Exp ',sprintf('%02i',projs(ex))]);
  for i=6:6
    fname=[pdir,sprintf('%02i',projs(ex)),'/his_',sprintf('%04i',i),'.nc'];
    ncid=netcdf.open(fname,'nowrite');
    for j=1:31
      %--- Read u and add to mean ---%
      u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 0 j-1],[401 482 36 1]);
      uu(1,:,:,:)=double(permute(u,[3 1 2]));
      meanu(ex,:,:,:)=meanu(ex,:,:,:)+uu; clear uu
 
      %--- Read v and add to mean ---%
      u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 0 j-1],[402 481 36 1]);
      uu(1,:,:,:)=double(permute(u,[3 1 2]));
      meanv(ex,:,:,:)=meanv(ex,:,:,:)+uu; clear uu
    end
    netcdf.close(ncid);
  end
end
%meanu=meanu./186; meanv=meanv./186;
%meanu=meanu./31;  meanv=meanv./31;
meanu(length(projs),:,:,:)=...
meanu(length(projs),:,:,:)./31;
meanv(length(projs),:,:,:)=...
meanv(length(projs),:,:,:)./31;
eval(['save ',datadir,'mean_JUNE_vels.mat meanu meanv -MAT'])

%load([datadir,'mean_JUNE_vels.mat']);

%=========================================================================%
%=                       Loop through to calculate EKE                   =%
%=========================================================================%
disp('Calculating EKE');

for ex=length(projs):length(projs)
  disp(['Exp ',sprintf('%02i',projs(ex))]);

  ke=zeros(36,402,482); count=0;
  for i=6:6
    fname=[pdir,sprintf('%02i',projs(ex)),'/his_',sprintf('%04i',i),'.nc'];
    ncid=netcdf.open(fname,'nowrite');
    for j=1:31

      %--- Read u and subtract mean ---%
      u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 0 j-1],[401 482 36 1]);
      u=double(permute(u,[3 1 2]))-squeeze(meanu(ex,:,:,:));
     
      %--- Read v and subtract mean ---%
      v=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 0 j-1],[402 481 36 1]);
      v=double(permute(v,[3 1 2]))-squeeze(meanv(ex,:,:,:));
      
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
    netcdf.close(ncid);
  end
  ke=ke./count;
  eval(['save ',datadir,'eke_JUNE_exp',sprintf('%02i',projs(ex)),...
        '.mat ke -MAT']);
end





