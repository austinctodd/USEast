%=============================================================================%
%=                                                                           =%
%=  Created by Austin C. Todd, NCSU, 2013                                    =%
%=                                                                           =%
%=  Program is designed to read in output from the US East model and calcu-  =%
%=  late the mean currents over several different seasons. Averages are      =%
%=  calculated for Spring (March - May), Summer (June-August), Fall          =%
%=  (September - November), and Winter (December - February).                =%
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
pdir='/gpfs_share/actodd/USeast-age/output/2005';

%=============================================================================%
%=                  Setup some grid and preliminary arrays                   =%
%=============================================================================%
refdate=datenum(2004,01,01); % Why the hell do we choose this date???

gridfile=[pdir,'/useast_his.nc'];

ncid=netcdf.open(gridfile,'nowrite');
h    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'         ));
lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'   ));
lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'   ));
mask1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_shelf'));
hmask=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'  ));
time =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'));
%netcdf.close(ncid)

lon=lon'; lat=lat';
nanmask1=mask1'./mask1';
nanmaskh=hmask'./hmask';
time=time./86400+refdate;

%=============================================================================%
%=        Loop to read in files and calculate surface mean velocities        =%
%=============================================================================%
refdate=datenum(2004,01,01); 

%--- Spring 2005 ---%
disp('Spring 2005')
inds=find(time>=datenum(2005,3,1) & time<datenum(2005,6,1));

count=0; mean_u=zeros(401,482); mean_v=zeros(402,481);
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 inds(i)-1],[401 482 1 1]);
  mean_u=mean_u+double(var);
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 inds(i)-1],[402 481 1 1]);
  mean_v=mean_v+double(var);
  count=count+1;
end
spring.u(1,:,:)=mean_u./count;
spring.v(1,:,:)=mean_v./count;

%--- Summer 2005 ---%
disp('Summer 2005');
inds=find(time>=datenum(2005,6,1) & time<datenum(2005,9,1));

count=0; mean_u=zeros(401,482); mean_v=zeros(402,481);
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 inds(i)-1],[401 482 1 1]);
  mean_u=mean_u+double(var);
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 inds(i)-1],[402 481 1 1]);
  mean_v=mean_v+double(var);
  count=count+1;
end
summ.u(1,:,:)=mean_u./count;
summ.v(1,:,:)=mean_v./count;

%--- Autumn 2005 ---%
disp('Autumn 2005');
inds=find(time>=datenum(2005,9 ,1) & time<datenum(2005,12,1));

count=0; mean_u=zeros(401,482); mean_v=zeros(402,481);
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 inds(i)-1],[401 482 1 1]);
  mean_u=mean_u+double(var);
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 inds(i)-1],[402 481 1 1]);
  mean_v=mean_v+double(var);
  count=count+1;
end
autumn.u(1,:,:)=mean_u./count;
autumn.v(1,:,:)=mean_v./count;

%--- Winter 2005/6 ---%
disp('Winter 2005/6');
inds=find(time>=datenum(2005,12,1) & time<datenum(2006, 3,1));

count=0; mean_u=zeros(401,482); mean_v=zeros(402,481);
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 inds(i)-1],[401 482 1 1]);
  mean_u=mean_u+double(var);
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 inds(i)-1],[402 481 1 1]);
  mean_v=mean_v+double(var);
  count=count+1;
end
winter.u(1,:,:)=mean_u./count;
winter.v(1,:,:)=mean_v./count;

%--- Spring 2006 ---%
disp('Spring 2006');
inds=find(time>=datenum(2006,3,1) & time<datenum(2006,6,1));

count=0; mean_u=zeros(401,482); mean_v=zeros(402,481);
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 inds(i)-1],[401 482 1 1]);
  mean_u=mean_u+double(var);
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 inds(i)-1],[402 481 1 1]);
  mean_v=mean_v+double(var);
  count=count+1;
end
spring.u(2,:,:)=mean_u./count;
spring.v(2,:,:)=mean_v./count;

%--- Summer 2006 ---%
disp('Summer 2006');
inds=find(time>=datenum(2006,6,1) & time<datenum(2006,9,1));

count=0; mean_u=zeros(401,482); mean_v=zeros(402,481);
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 inds(i)-1],[401 482 1 1]);
  mean_u=mean_u+double(var);
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 inds(i)-1],[402 481 1 1]);
  mean_v=mean_v+double(var);
  count=count+1;
end
summ.u(2,:,:)=mean_u./count;
summ.v(2,:,:)=mean_v./count;

%--- Autumn 2006 ---%
inds=find(time>=datenum(2006,9 ,1) & time<datenum(2006,12,1));

count=0; mean_u=zeros(401,482); mean_v=zeros(402,481);
for i=1:length(inds)
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 inds(i)-1],[401 482 1 1]);
  mean_u=mean_u+double(var);
  var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 inds(i)-1],[402 481 1 1]);
  mean_v=mean_v+double(var);
  count=count+1;
end
autumn.u(2,:,:)=mean_u./count;
autumn.v(2,:,:)=mean_v./count;

clear mean_u mean_v count;

netcdf.close(ncid);

save meanvels.mat spring summ autumn winter -MAT


return;





%--- Spring 2005 ---%
stind=find(time==datenum(2005,3,1))-1;
enind=find(time==datenum(2005,6,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 stind],[401 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 stind],[402 481 1 length(stind:enind-1)]);

spring.mu(1,:,:)=squeeze(nanmean(var1,3));
spring.mv(1,:,:)=squeeze(nanmean(var2,3));
clear var1 var2 ;

%--- Summer 2005 ---%
stind=find(time==datenum(2005,6,1))-1;
enind=find(time==datenum(2005,9,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 stind],[401 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 stind],[402 481 1 length(stind:enind-1)]);
summ.mu(1,:,:)=squeeze(nanmean(u,3));
summ.mv(1,:,:)=squeeze(nanmean(v,3));

clear var1 var2 ;

%--- Autumn 2005 ---%
stind=find(time==datenum(2005,9 ,1))-1;
enind=find(time==datenum(2005,12,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 stind],[401 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 stind],[402 481 1 length(stind:enind-1)]);
autumn.mu(1,:,:)=squeeze(nanmean(u,3));
autumn.mv(1,:,:)=squeeze(nanmean(v,3));

clear var1 var2;

%--- Winter 2005/6 ---%
stind=find(time==datenum(2005,12,1))-1;
enind=find(time==datenum(2006, 3,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 stind],[401 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 stind],[402 481 1 length(stind:enind-1)]);
winter.mu(1,:,:)=squeeze(nanmean(u,3));
winter.mv(1,:,:)=squeeze(nanmean(v,3));

clear var1 var2;

%--- Spring 2006 ---%
stind=find(time==datenum(2006,3,1))-1;
enind=find(time==datenum(2006,6,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 stind],[401 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 stind],[402 481 1 length(stind:enind-1)]);
spring.mu(2,:,:)=squeeze(nanmean(u,3));
spring.mv(2,:,:)=squeeze(nanmean(v,3));

clear var1 var2;

%--- Summer 2006 ---%
stind=find(time==datenum(2006,6,1))-1;
enind=find(time==datenum(2006,9,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 stind],[401 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 stind],[402 481 1 length(stind:enind-1)]);
summ.mu(2,:,:)=squeeze(nanmean(u,3));
summ.mv(2,:,:)=squeeze(nanmean(v,3));

clear var1 var2;

%--- Autumn 2006 ---%
stind=find(time==datenum(2006,9 ,1))-1;
enind=find(time==datenum(2006,12,1))-1;

var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),...
                  [0 0 35 stind],[401 482 1 length(stind:enind-1)]);
var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                  [0 0 35 stind],[402 481 1 length(stind:enind-1)]);
autumn.mu(2,:,:)=squeeze(nanmean(u,3));
autumn.mv(2,:,:)=squeeze(nanmean(v,3));

clear var1 var2;



