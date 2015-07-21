%======================== plot_tide_comparison.m =========================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  Program obtains all the same information as plot_tides.m, except the =%
%=  ROMS output is on omgsrv2. So, the data for CFSR and for the tide    =%
%=  gauges is just saved into a MAT file to plot there.                  =%
%=========================================================================%

%=========================================================================%
%=                        Load dependent libraries                       =%
%=========================================================================%
addpath(genpath('/he_data/he/zhigang/Matlab_Codes/mexcdf2009/'));
addpath('/home/actodd/MYMATLAB/');
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));

%=========================================================================%
%=                         User-defined indexes                          =%
%=========================================================================%
tStart = datenum(2004,01,01);
tEnd   = datenum(2006,01,01);

%=========================================================================%
%=                   Set ROMS file name and load data                    =%
%=========================================================================%
roms.dir ='/gpfs_share/rhe/actodd/USeast-age/output/OBC/2004_2/';
roms.f   =[roms.dir,'useast_his.nc'];
roms.lon =nc_varget(roms.f,'lon_rho'                 );
roms.lat =nc_varget(roms.f,'lat_rho'                 );
roms.mask=nc_varget(roms.f,'mask_rho'                );

%=========================================================================%
%=                   Set CFSR file name and load data                    =%
%=========================================================================%
cfsr.f='/he_data/he/actodd/ROMS-age/Data/frc/USeast-frc.nc';
cfsr.time=nc_varget(cfsr.f,'time',1443,736,3);
cfsr.date=cfsr.time+datenum(2004,1,1);
cfsr.lon =nc_varget(cfsr.f,'lon');
cfsr.lat =nc_varget(cfsr.f,'lat');

%=========================================================================%
%=                      Set list of all station IDs                      =%
%=========================================================================%
List4sID=[8410140;8413320;8418150;8443970;8449130;8447930;8452660
          8461490;8465705;8516945;8510560;8531680;8534720;8557380
          8570283;8632200;8575512;8636580;8651370;8656483;8661070
          8665530;8670870;8720218;8721604;8723214;8724580;8723970
          8725110;8726384;8727520;8728690;8729108;8735180;8760922
          8770570;8771341;8775870;2695540];
nSta = length(List4sID);

% Tide station information matrix
load('Data/US_TideInfo.mat');

%=========================================================================%
%=                   Set Tide file name and load data                    =%
%=========================================================================%
for i=1:nSta
  tf=['/he_data/he/actodd/DATA/tides/',sprintf('%7i',List4sID(i)),...
          '.dat'];
  tdata=load(tf);
  tide.date{i}=datenum(tdata(:,1),tdata(:,2),tdata(:,3),...
                       tdata(:,4),tdata(:,5),tdata(:,5).*0);
  tide.data{i}=tdata(:,6);
  tide.id{i}  =List4sID(i);
  
  [rxind(i),ryind(i)]=find_closest_roms_point(roms,ss,List4sID(i));  
  [xind,yind]=find_closest_cfsr_point(cfsr,ss,List4sID(i));
  cfsr.pres(i,:)=nc_varget(cfsr.f,'Pair',[1443 yind xind],[736 1 1],[3 1 1]); 
    
end

save('data_for_srv2.mat','tide','rxind','ryind','cfsr','-MAT');


