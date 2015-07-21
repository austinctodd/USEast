%=========================== make_USeast_ads.m ===========================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  Program creates Adjoint Sensitivity and Scope files for the US East  =%
%=  Coast water residence time model. Values are 1 inside control volume =%
%=  and 0 elsewhere.                                                     =%
%=                                                                       =%
%=========================================================================%

%=========================================================================%
%=                        Load dependent libraries                       =%
%=========================================================================%
addpath(genpath('/he_data/he/zhigang/Matlab_Codes/mexcdf2009/'));
addpath('/home/actodd/MYMATLAB/');
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));

%=========================================================================%
%=                Set filenames for ADS and Scope files                  =%
%=========================================================================%
adsfile  ='/he_data/he/actodd/ROMS-age/Data/rtime/USeast-ads.nc';
scopefile='/he_data/he/actodd/ROMS-age/Data/rtime/USeast-scope.nc';

%=========================================================================%
%=                  Read US East ROMS grid infomation                    =%
%=========================================================================%
roms.f   ='/he_data/he/actodd/ROMS-age/Data/grid/USeast-grid-fundy.nc';
roms.lon =nc_varget(roms.f,'lon_rho'                 );
roms.lat =nc_varget(roms.f,'lat_rho'                 );
roms.mask=nc_varget(roms.f,'mask_shelf'              );

%=========================================================================%
%=       Make u and v scopes and set ads values at all depths/times      =%
%=========================================================================%
[roms.masku,roms.maskv,roms.maskp]=uvp_masks(roms.mask');
ads=zeros(2,36,482,402);
for i=1:36
  ads(1,i,:,:)=roms.mask;
  ads(2,i,:,:)=roms.mask;
end

%=========================================================================%
%=                Write data to Adjoint Sensitivity File                 =%
%=========================================================================%
nc_varput(adsfile,'ocean_time',[0 100000000],0,2);
nc_varput(adsfile,'age_01',ads);
nc_varput(adsfile,'age_02',ads);
nc_varput(adsfile,'age_03',ads);
nc_varput(adsfile,'age_04',ads);
nc_varput(adsfile,'age_05',ads);
nc_varput(adsfile,'age_06',ads);

%=========================================================================%
%=                     Write data to Grid Scope File                     =%
%=========================================================================%
nc_varput(scopefile,'scope_rho',roms.mask );
nc_varput(scopefile,'scope_u',  roms.masku');
nc_varput(scopefile,'scope_v',  roms.maskv');

