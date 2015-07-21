%=============================================================================%
%=                                                                           =%
%=  Created by Austin C. Todd, NCSU, 2013                                    =%
%=                                                                           =%
%=  Program is designed to read in age data from the US East model output    =%
%=  and calculate the mean age over several different seasons. The mean age  =%
%=  will then be averaged over several time periods to correspond with calc- =%
%=  ulations done in Zhang et al., 2010 near the Hudson Bay shelf water      =%
%=  region. Averages will be for Spring (March - May), Summer (June-August), =%
%=  Fall (September - November), and Winter (December - February).           =%
%=                                                                           =%
%=============================================================================%

%--- Startup and setting of the path ---%
addpath('~/MYMATLAB/');
addpath(genpath('~/MYMATLAB/ROMS-matlab/'));
addpath(genpath('/Users/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=============================================================================%
%=                                                                           =%
%=       !!!!!!!!!!!        USER-DEFINED VARIABLES      !!!!!!!!!!           =%
%=                                                                           =%
%=============================================================================%
tracer='shelf'; % 'shelf' or 'deep' or 'river'
pdir='/Volumes/Black_box/Data/USeast-age/output/clim/averages/';

%=============================================================================%
%=                  Setup some grid and preliminary arrays                   =%
%=============================================================================%
refdate=datenum(2004,01,01); % Why the hell do we choose this date???

gridfile=['/Volumes/Black_box/Data/USeast-age/output/rtime/useast_his.nc'];

ncid=netcdf.open(gridfile,'nowrite');
h    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'         ));
lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'   ));
lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'   ));
mask1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_shelf'));
hmask=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'  ));
%time =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'));
netcdf.close(ncid)

nc=netcdf.open([pdir,'avg_3hrly.nc'],'nowrite');

lon=lon'; lat=lat';
nanmask1=mask1'./mask1';
nanmaskh=hmask'./hmask';
%time=time./86400+refdate;

%=============================================================================%
%=           Loop to read in files and calculate surface mean age            =%
%=============================================================================%
refdate=datenum(2004,01,01); 

day_count=[31,28,31,30,31,30,31,31,30,31,30,31];

seasonal_u   =zeros(401,482,35,4);
seasonal_v   =zeros(402,481,35,4);
seasonal_age1=zeros(402,482,35,4);
seasonal_age2=zeros(402,482,35,4);
seasonal_age3=zeros(402,482,35,4);

for season=1:4
  disp(['Season ',sprintf('%1i',season)])
  if (season==1)
    stind=0;
  else
    stind=(sum(day_count(1:(season-1)*3))-1)*8;
  end
  count=sum(day_count(((season-1)*3+1):((season-1)*3+3)))*8;

  for i=1:35
    fprintf(1,['Depth ',sprintf('%02i',i),': ']);
 
    %--- Get u velocities ---%
    fprintf(1,'U...');
    u=netcdf.getVar(nc,netcdf.inqVarID(nc,'u'),...
                    [0 0 i-1 stind],[401 482 1 count]);
    seasonal_u(:,:,i,season)=nanmean(double(squeeze(u)),3); 
    clear u;
  
    %--- Get v velocities ---%
    fprintf(1,'V...');
    u=netcdf.getVar(nc,netcdf.inqVarID(nc,'v'),...
                      [0 0 i-1 stind],[402 481 1 count]);
    seasonal_v(:,:,i,season)=nanmean(double(squeeze(u)),3);
    clear u
  
    %--- Deep ocean source water ---%
    fprintf(1,'Age 1...');
    meanage=netcdf.getVar(nc,netcdf.inqVarID(nc,'mean_age_01'),...
                      [0 0 i-1 stind],[402 482 1 count]);
    seasonal_age1(:,:,i,season) = nanmean(double(squeeze(meanage)),3);
    clear meanage
  
    %--- Mississippi River source water ---%
    fprintf(1,'Age 2...');
    meanage=netcdf.getVar(nc,netcdf.inqVarID(nc,'mean_age_02'),...
                      [0 0 i-1 stind],[402 482 1 count]);
    seasonal_age2(:,:,i,season) = nanmean(double(squeeze(meanage)),3);       
    clear meanage

    %--- All river source water ---%
    fprintf(1,'Age 3...');
    meanage=netcdf.getVar(nc,netcdf.inqVarID(nc,'mean_age_03'),...
                      [0 0 i-1 stind],[402 482 1 count]);
    seasonal_age3(:,:,i,season) = nanmean(double(squeeze(meanage)),3);       
    clear meanage

    fprintf(1,'\n');
  end
end

nanmask=nanmaskh;
netcdf.close(nc);
return;
