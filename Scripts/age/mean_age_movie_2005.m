%=============================================================================%
%=                                                                           =%
%=  Created by Austin C. Todd, NCSU, 2013                                    =%
%=                                                                           =%
%=  Program is designed to read in age information from the US East model    =%
%=  and calculate the mean age.  The model has three different "source"      =%
%=  water regions:  1. Shelf water, 2. River water, 3. Deep ocean water.     =%
%=                                                                           =%
%=  For each source region, there are two tracer types:                      =%
%=             C = water age tracer concentration                            =%
%=         alpha = age concentration                                         =%
%=                                                                           =%
%=  where                                                                    =%
%=                                                                           =%
%=         dC/dt = P - D - del*(UC       - K*(delC      ))                   =%
%=   d(alpha)/dt =     C - del*(U(alpha) - K*(del(alpha)))                   =%
%=                                                                           =%
%=  P = source, D = sink, del = d/dx+d/dy, * = dot product, U = velocity,    =%
%=  and K = diffusivity                                                      =%
%=                                                                           =%
%=                      The mean age A = alpha / C                           =%
%=                                                                           =%
%=============================================================================%

%--- Startup and setting of the path ---%
addpath('~/MYMATLAB/');
addpath(genpath('~/MYMATLAB/ROMS-matlab/'));
addpath(genpath('/raid0/xzeng2/software/mexcdf/'))
format long g; format compact;

%=============================================================================%
%=                                                                           =%
%=       !!!!!!!!!!!        USER-DEFINED VARIABLES      !!!!!!!!!!           =%
%=                                                                           =%
%=============================================================================%
tracer='shelf'; % 'shelf' or 'deep' or 'river'
pdir='/raid0/actodd/ROMS-exp/USeast-age/output/2005';

%=============================================================================%
%=                  Setup some grid and preliminary arrays                   =%
%=============================================================================%
refdate=datenum(2004,01,01); % Why the hell do we choose this date???

gridfile=[pdir,'/USeast-his_0001.nc'];
%gridfile='/raid0/actodd/ROMS-exp/USeast-age/output/test3_2month/USeast_his_0060.nc';
%Vstretching=2;  Vtransform=2;
%theta_s=nc_varget(gridfile,'theta_s'   );
%theta_b=nc_varget(gridfile,'theta_b'   );
%Tcline =nc_varget(gridfile,'Tcline'    );
%hc     =nc_varget(gridfile,'hc'        );
ncid=netcdf.open(gridfile,'nowrite');
h    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'         ));
lon  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'   ));
lat  =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'   ));
mask1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_shelf'));
hmask=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'  ));
netcdf.close(ncid)
%z      =set_depth(2,2,theta_s,theta_b,Tcline,36,1,h,h.*0,0);
%zw     =set_depth(2,2,theta_s,theta_b,Tcline,36,5,h,h.*0,0);

%z=-1.*permute(z,[3 1 2]); zw=-1.*permute(zw,[3 1 2]);
%Hz=diff(zw,1);

%for i=1:36
% zmask1(i,:,:)=mask1; zmask2(i,:,:)=mask2;, zmask3(i,:,:)=mask3;
%end
%zmiss=z; zmiss(find(zmiss>200 ))=nan;
%zmiss=zmiss.*0+1;

%zdeep=z; zdeep( find(zdeep <2500))=nan;
%zdeep=zdeep.*0+1;

lon=lon'; lat=lat';
nanmask1=mask1'./mask1';
nanmaskh=hmask'./hmask';

%=============================================================================%
%=        Load in colorbar and make one that is just shades of purple        =%
%=============================================================================%
cmap=load('/raid0/actodd/LTRANS-utils/density_colormap.txt');
cmap(end,:)=cmap(end-1,:);
for i=29:-1:1
  cmap(i,1)=cmap(i+1,1)+(0.9-cmap(30,1))/30;
  cmap(i,2)=cmap(i+1,2)+(0.9/30);
  cmap(i,3)=cmap(i+1,3)+(0.9-cmap(30,3))/30;
end
cmap=cmap(end:-1:1,:);
%cmap(end,:)=[1 1 1];

%=============================================================================%
%=           Loop to read in files and plot surface age data                 =%
%=============================================================================%
%refdate=datenum(1858,11,17); % Why the hell do we choose this date???
refdate=datenum(2004,01,01); % Why the hell do we choose this date???

count=1;
%fname=[pdir,'_mpdata/USeast_daily.nc'];
for day=0:200

  if (day<=30) 
    fnum=day; 
    fname=[pdir,'/USeast-his_0001.nc'];
  else 
    fnum=mod(day-1,30);
    fname=[pdir,'/USeast-his_',sprintf('%04i',floor((day-1)/30)+1),'.nc'];
  end

  ncid=netcdf.open(fname,'nowrite');

  ocean_time=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'),...
                           [fnum],[1]);

  disp([datestr(datenum(2004,01,1+ocean_time/86400))]);

  if strcmp(tracer,'deep')
    var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_01'),...
                      [0 0 35 fnum],[402 482 1 1]);
    var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_02'),...
                      [0 0 35 fnum],[402 482 1 1]);
    nanmask=nanmask1;
  elseif strcmp(tracer,'Miss');
    var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_03'),...
                      [0 0 35 fnum],[402 482 1 1]);
    var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_04'),...
                      [0 0 35 fnum],[402 482 1 1]);
    nanmask=nanmaskh;
  else
    var1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_05'),...
                      [0 0 35 fnum],[402 482 1 1]);
    var2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_06'),...
                      [0 0 35 fnum],[402 482 1 1]);
    nanmask=nanmaskh;
  end
  netcdf.close(ncid);
  var1=double(var1); var2=double(var2);

%  var1(find(var1 < 1e-5))=0; 
var1=var1';
%  var2(find(var1 < 1e-5))=0; 
var2=var2';
  mean_age=var2./var1;
  mean_age(find(mean_age < 0 ))=nan;
  mean_age(find(mean_age > day))=day;  
 
  %--- PLOT ---%
  figure(1); clf;
  pcolor(lon,lat,mean_age.*nanmask); shading flat; axis xy; axis image;
  hold on;
  contour(lon,lat,hmask',[0 0],'k');
  caxis([0 60]);
  axis([-98 -80 18 31]);
  set(gca,'ytick',[20:2:30]);

  %--- Make Title for entire figure ---%
  hTitle=title(['Mean age (days) :: Model integration day ',...
               sprintf('%02i',day)],'fontsize',14,'fontweight','bold');

  colormap(1-pink);
%  colormap(jet)

  %--- Include Colorbar ---%
  hCbar=colorbar('southoutside');%,'position',[0.05 0.24 0.43 0.05]);

  drawnow;
  wysiwyg;
  pause(0.25)

%  print('-depsc2','-painters',['/raid0/actodd/PLOTS/USeast-age/frames/',...
%                               'meanage_',sprintf('%02i',day),'.eps']);

end


