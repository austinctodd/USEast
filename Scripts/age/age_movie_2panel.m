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
%                !!!!!!!!    USER-DEFINED VARIABLES    !!!!!!!!              =%
%=                                                                           =%
%=============================================================================%
tracer='shelf'; %choose either 'shelf', 'deep', or 'river'
pdir ='/raid0/actodd/ROMS-exp/USeast-age/output/test3_2month';

%=============================================================================%
%=                  Setup some grid and preliminary arrays                   =%
%=============================================================================%
refdate=datenum(1858,11,17); % Why the hell do we choose this date???

gridfile=[pdir,'/USeast_his_0060.nc'];
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
hmask  =nc_varget(gridfile,'mask_rho'  );
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

if strcmp(tracer,'shelf')
  nanmask2=mask1./mask1;
elseif strcmp(tracer,'deep')
  nanmask2=mask3./mask3;
else
  nanmask2=mask2/mask2;
end
nanmask =hmask./hmask;

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
%cmap(1,:)=[1 1 1];

cmap2(:,1)=interp1([0:29],cmap(1:30,1),[0:29/255:29]);
cmap2(:,2)=interp1([0:29],cmap(1:30,2),[0:29/255:29]);
cmap2(:,3)=interp1([0:29],cmap(1:30,3),[0:29/255:29]);
cmap2(1,:)=[1 1 1];

%=============================================================================%
%=           Loop to read in files and plot surface age data                 =%
%=============================================================================%
refdate=datenum(1858,11,17); % Why the hell do we choose this date???

count=1;
%pdir='/raid0/actodd/ROMS-exp/USeast-age/output/test3_41days/';
fname=[pdir,'_mpdata/USeast_daily.nc'];
for day=1:61
  ocean_time=nc_varget(fname,'ocean_time');
  disp(['Day : ',sprintf('%02i',day)]);

  if (strcmp(tracer,'shelf'))
    var1=nc_varget(fname,['age_01'],[day-1 35 0 0],[1 1 482 402]);
    var2=nc_varget(fname,['age_02'],[day-1 35 0 0],[1 1 482 402]);
  elseif (strcmp(tracer,'deep'))
    var1=nc_varget(fname,['age_05'],[day-1 0  0 0],[1 1 482 402]);
    var2=nc_varget(fname,['age_06'],[day-1 0  0 0],[1 1 482 402]);
  else
    var1=nc_varget(fname,['age_03'],[day-1 35 0 0],[1 1 482 402]);
    var2=nc_varget(fname,['age_04'],[day-1 35 0 0],[1 1 482 402]);
  end
  mean_age=var2./var1;
  
  %--- PLOT ---%
  figure(1); clf;
  ax1=axes('position',[0.04 0.15 0.45 0.8]);
  pcolor(lon,lat,var1.*nanmask); shading flat; axis xy; axis image;
  hold on;
  contour(lon,lat,hmask,[0 0],'k');
  caxis([0 1]);
  axis([-98 -80 18 31]);
  set(gca,'ytick',[20:2:30]);

  ax2=axes('position',[0.52 0.15 0.45 0.8]);
  pcolor(lon,lat,var2.*nanmask2); shading flat; axis xy; axis image;
  hold on;
  contour(lon,lat,hmask,[0 0],'k');
  caxis([0 30])
  axis([-98 -80 18 31]);
  set(gca,'ytick',[20:2:30]);  

  %--- Make Title for entire figure ---%
  ax3=axes('position',[0.05 0.75 0.9 0.1]);
  hText=text(0.5,0.5,['Day ',sprintf('%02i',day)],'fontsize',14,...
             'fontweight','bold','horizontalalignment','center');
  axis off;

  %--- Label tracer concentration panel ---%
  ax_Ctitle=axes('position',[0.3 0.334 0.19 0.1]);
  axis off;
  hCTitle1=text(.98,.45,      'tracer','HorizontalAlignment','right');
  hCTitle2=text(.98,.2,'concentration','HorizontalAlignment','right');
  set([hCTitle1,hCTitle2],'fontsize',11);

  %--- Label tracer concentration panel ---%
  ax_Atitle=axes('position',[0.78 0.334 0.19 0.1]);
  axis off;
  hATitle1=text(.98,.45,   'water age','HorizontalAlignment','right');
  hATitle2=text(.98,.2,'concentration','HorizontalAlignment','right');
  set([hATitle1,hATitle2],'fontsize',11);

  colormap(cmap2);

  %===========================================================================%
  %                            Make tracer colorbar                           %
  %===========================================================================%
  ax_CCbar=axes('position',[0.05 0.24 0.43 0.05]);
  xbar=[0      0:1/67:1        1    1:-1/67:0      ];
  ybar=[0     [0:1/67:1].*0+1  1   [1:-1/67:0].*0+0];
  cbar=[0      0:1/67:1        1    1:-1/67:0      ];
  xbarlines=[0 .2 .4 .6 .8 1];
  xbarticks={'0';'0.2';'0.4';'0.6';'0.8';'1'};

  %--- Make our own colorbar in the land-masked portion of the domain ---%
  hCbar=patch(xbar,ybar,cbar); hold on; axis off;
  set(hCbar,'linewidth',1);

  %--- Draw colorbar tick lines ---%
  for j=2:length(xbarlines)-1
    plot([xbarlines(j) xbarlines(j)],[0 0.2],'k','linewidth',1);
    plot([xbarlines(j) xbarlines(j)],[1 0.8],'k','linewidth',1);
  end

  %--- Write values at colorbar tick lines ---%
  for j=1:length(xbarlines)
     hTicks{j}=text(xbarlines(j),-0.3,xbarticks(j),...
                    'FontSize',10,'HorizontalAlignment','Center');
  end

  %===========================================================================%
  %                            Make tracer colorbar                           %
  %===========================================================================%
  ax_ACbar=axes('position',[0.53 0.24 0.43 0.05]);
  xbarlines=[0 .2 .4 .6 .8 1];
  xbarticks={'0';'6';'12';'18';'24';'30'};

  %--- Make our own colorbar in the land-masked portion of the domain ---%
  hCbar=patch(xbar,ybar,cbar); hold on; axis off;
  set(hCbar,'linewidth',1);

  %--- Draw colorbar tick lines ---%
  for j=2:length(xbarlines)-1
    plot([xbarlines(j) xbarlines(j)],[0 0.2],'k','linewidth',1);
    plot([xbarlines(j) xbarlines(j)],[1 0.8],'k','linewidth',1);
  end

  %--- Write values at colorbar tick lines ---%
  for j=1:length(xbarlines)
     hTicks{j}=text(xbarlines(j),-0.3,xbarticks(j),...
                    'FontSize',10,'HorizontalAlignment','Center');
  end

  drawnow;
  wysiwyg;
  pause(0.25)
%return;
 % print('-depsc2','-painters',['/raid0/actodd/PLOTS/USeast-age/frames/',...
 %                              'day_',sprintf('%02i',day),'.eps']);

end


