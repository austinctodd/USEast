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
%=                  Setup some grid and preliminary arrays                   =%
%=============================================================================%
refdate=datenum(1858,11,17); % Why the hell do we choose this date???
pdir='/raid0/actodd/ROMS-exp/USeast-age/output/test3_2month/';

gridfile=[pdir,'USeast_his_0060.nc'];
Vstretching=2;  Vtransform=2;
theta_s=nc_varget(gridfile,'theta_s'   );
theta_b=nc_varget(gridfile,'theta_b'   );
Tcline =nc_varget(gridfile,'Tcline'    );
hc     =nc_varget(gridfile,'hc'        );
h      =nc_varget(gridfile,'h'         );
lon    =nc_varget(gridfile,'lon_rho'   );
lat    =nc_varget(gridfile,'lat_rho'   );
mask1  =nc_varget(gridfile,'mask_shelf');
mask2  =nc_varget(gridfile,'mask_Miss' );
mask3  =nc_varget(gridfile,'mask_deep' );
z      =set_depth(2,2,theta_s,theta_b,Tcline,36,1,h,h.*0,0);
zw     =set_depth(2,2,theta_s,theta_b,Tcline,36,5,h,h.*0,0);

z=-1.*permute(z,[3 1 2]); zw=-1.*permute(zw,[3 1 2]);
Hz=diff(zw,1);

for i=1:36
 zmask1(i,:,:)=mask1; zmask2(i,:,:)=mask2;, zmask3(i,:,:)=mask3;
end
zmiss=z; zmiss(find(zmiss>200 ))=nan;
zmiss=zmiss.*0+1;

zdeep=z; zdeep( find(zdeep <2500))=nan;
zdeep=zdeep.*0+1;

%=============================================================================%
%=                        Loop to read in files                              =%
%=============================================================================%
refdate=datenum(1858,11,17); % Why the hell do we choose this date???

count=1;
!for day=1:30
day=60;
  fname=[pdir,'USeast_his_',sprintf('%04i',day),'.nc'];
  ocean_time=nc_varget(fname,'ocean_time');
  disp(['Day : ',sprintf('%02i',day)]);
!  for i=1:length(ocean_time);
i=8;
    for age=1:3
      var=nc_varget(fname,['age_',sprintf('%02i',age*2-1)],[i-1 0 0 0],[1 36 482 402]);
      if (age==2) var=var.*zmiss; elseif (age==3) var=var.*zdeep; end
      eval(['var1=squeeze(nansum(var.*Hz.*zmask',sprintf('%1i',ceil(age)),...
                                       './zmask',sprintf('%1i',ceil(age)),'))./h;']);
      var=nc_varget(fname,['age_',sprintf('%02i',age*2  )],[i-1 0 0 0],[1 36 482 402]);
      if (age==2) var=var.*zmiss; elseif (age==3) var=var.*zdeep; end
      eval(['var2=squeeze(nansum(var.*Hz.*zmask',sprintf('%1i',ceil(age)),...
                                       './zmask',sprintf('%1i',ceil(age)),'))./h;']);
!      eval(['mean_age_',sprintf('%1i',age),'(count,:,:)=var2./var1;']);
      eval(['mean_age_',sprintf('%1i',age),'=var2./var1;']);
!    end
    count=count+1;
  end
!end

