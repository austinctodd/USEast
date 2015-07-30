%=============================================================================%
%=                                                                           =%
%=  Created by Austin C. Todd, NCSU, 2013                                    =%
%=                                                                           =%
%=  Program is designed to read residence time (RT) from the US East model   =%
%=  and calculate the mean RT over several different seasons. The mean RT    =%
%=  will then be averaged over several time periods to correspond with calc- =%
%=  ulations done in Zhang et al., 2010 near the Hudson Bay shelf water      =%
%=  region. Averages will be for Spring (March - May), Summer (June-August), =%
%=  Fall (September - November), and Winter (December - February).           =%
%=                                                                           =%
%=============================================================================%

%--- Startup and setting of the path ---%
addpath('~/MYMATLAB/');
addpath(genpath('~/MYMATLAB/ROMS-matlab/'));
addpath(genpath('/home/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=============================================================================%
%=                                                                           =%
%=       !!!!!!!!!!!        USER-DEFINED VARIABLES      !!!!!!!!!!           =%
%=                                                                           =%
%=============================================================================%
tracer='shelf'; % 'shelf' or 'deep' or 'river'
pdir='/gpfs_share/actodd/USeast-rtime/output/compressed/';

%=============================================================================%
%=                  Setup some grid and preliminary arrays                   =%
%=============================================================================%
refdate=datenum(2004,01,01); % Why the hell do we choose this date???

gridfile=['/gpfs_share/actodd/USeast-age/output/2005/useast_his.nc'];

ncid=netcdf.open(gridfile,'nowrite');
h    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'         ));
lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'   ));
lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'   ));
mask1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_shelf'));
hmask=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'  ));
netcdf.close(ncid)

lon=lon'; lat=lat';
nanmask1=mask1'./mask1';
nanmaskh=hmask'./hmask';

ncid=netcdf.open([pdir,'useast_adj.nc'],'nowrite');
time =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'));
time=time./86400+refdate;

%=============================================================================%
%=      Loop to read in files and calculate surface mean residence time      =%
%=============================================================================%
refdate=datenum(2004,01,01); 

%--- Spring 2005 ---%
disp('Spring 2005')
inds=find(time>=datenum(2005,3,1) & time<datenum(2005,6,1));

count=0; mean_rt=double(hmask).*0;
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 inds(i)-1],[402 482 1 1]);
  mean_rt=mean_rt+double(var)./86400;
  count=count+1;
end
spring.rt(1,:,:)=mean_rt./count;

%--- Summer 2005 ---%
disp('Summer 2005');
inds=find(time>=datenum(2005,6,1) & time<datenum(2005,9,1));

count=0; mean_rt=double(hmask).*0;
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 inds(i)-1],[402 482 1 1]);
  mean_rt=mean_rt+double(var)./86400;
  count=count+1;
end
summ.rt(1,:,:)=mean_rt./count;

%--- Autumn 2005 ---%
disp('Autumn 2005');
inds=find(time>=datenum(2005,9 ,1) & time<datenum(2005,12,1));

count=0; mean_rt=double(hmask).*0;
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 inds(i)-1],[402 482 1 1]);
  mean_rt=mean_rt+double(var)./86400;
  count=count+1;
end
autumn.rt(1,:,:)=mean_rt./count;

%--- Winter 2005/6 ---%
disp('Winter 2005/6');
inds=find(time>=datenum(2005,12,1) & time<datenum(2006, 3,1));

count=0; mean_rt=double(hmask).*0;
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 inds(i)-1],[402 482 1 1]);
  mean_rt=mean_rt+double(var)./86400;
  count=count+1;
end
winter.rt(1,:,:)=mean_rt./count;

%--- Spring 2006 ---%
disp('Spring 2006');
inds=find(time>=datenum(2006,3,1) & time<datenum(2006,6,1));

count=0; mean_rt=double(hmask).*0;
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 inds(i)-1],[402 482 1 1]);
  mean_rt=mean_rt+double(var)./86400;
  count=count+1;
end
spring.rt(2,:,:)=mean_rt./count;

%--- Summer 2006 ---%
disp('Summer 2006');
inds=find(time>=datenum(2006,6,1) & time<datenum(2006,9,1));

count=0; mean_rt=double(hmask).*0;
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 inds(i)-1],[402 482 1 1]);
  mean_rt=mean_rt+double(var)./86400;
  count=count+1;
end
summ.rt(2,:,:)=mean_rt./count;

%--- Autumn 2006 ---%
inds=find(time>=datenum(2006,9 ,1) & time<datenum(2006,12,1));

count=0; mean_rt=double(hmask).*0;
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                  [0 0 35 inds(i)-1],[402 482 1 1]);
  mean_rt=mean_rt+double(var)./86400;
  count=count+1;
end
autumn.rt(2,:,:)=mean_rt./count;

netcdf.close(ncid);

save meanrt.mat spring summ autumn winter -MAT

