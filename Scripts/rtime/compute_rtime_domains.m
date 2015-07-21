%=============================================================================%
%=                                                                           =%
%=  Created by Austin C. Todd, NCSU, 2015                                    =%
%=                                                                           =%
%=  Program reads in output from the US East residence time model to compute =%
%=  the residence time within each various sub-region.  The residence time   =%
%=  is a function of the volume of the control domain, so must be split      =%
%=  independently for each of the four different domains (GOM,MAB,SAB,NIS.   =%
%=                                                                           =%
%=============================================================================%

%--- Startup and setting of the path ---%
addpath('~/MYMATLAB/');
addpath(genpath('~/MYMATLAB/ROMS-matlab/'));
addpath(genpath('/Users/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=============================================================================%
%=                  Setup some grid and preliminary arrays                   =%
%=============================================================================%
gridfile=['/Volumes/Black_box/Data/USeast/Data/grd/USeast-grid.nc'];

ncid=netcdf.open(gridfile,'nowrite'); 
  h          =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'          ));
  lon        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'    ));
  lat        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'    ));
  mask       =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_shelf' ));
  hmask      =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'   ));
  scope      =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'scope_rho'  ));
  pm         =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'pm'         ));
  pn         =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'pn'         ));  
netcdf.close(ncid)

%--- Various Grid Parameters ---%
Vtransform = 2 ; Vstretching = 2 ;
theta_s = 7 ; theta_b = 0.1 ;
hc = 200 ; N=36;

%=============================================================================%
%=                     Read in grid adjoint information                      =%
%=============================================================================%
adsfile=['/Volumes/Black_box/Data/USeast/Data/grd/USeast-ads.nc' ];
ncid=netcdf.open(adsfile,'nowrite'); 
  funal      =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_02'     ),...
                            [0 0 0 0],[402 482 36 1]);
  funal = permute(funal,[3 1 2]);
netcdf.close(ncid)
return;
%=============================================================================%
%=                     Calculate various grid variables                      =%
%=============================================================================%
z_w = set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N,5,h,h.*0,0);
z_w = permute(z_w,[3 1 2]);
dxdy     = 1./pm./pn;
dxdy     = reshape(dxdy,[1 size(dxdy)]);
dxdy     = repmat(dxdy,[N 1]);
dz       = z_w(2:end,:,:) - z_w(1:end-1,:,:);
vol      = dz .* dxdy;
scope    = scope';
scope3d  = reshape(scope,[1 size(scope)]);
scope3d  = repmat(scope3d,[N 1]);
scope3d  = scope3d .* funal;
nvol     = vol .* scope3d;
ttvol    = sum(sum(sum(nvol)));
numgrd   = length(find(scope3d > 0.5));
meanvol  = ttvol ./ numgrd;
grdfct   = meanvol ./ vol;

%--- SAB ---%
sab.xs = 182; sab.xe = 243; 
sab.ys = 205; sab.ye = 315;
sab.vol    =vol(    :,sab.xs:sab.xe,sab.ys:sab.ye);
sab.funal  =funal(  :,sab.xs:sab.xe,sab.ys:sab.ye);
sab.scope3d=scope3d(:,sab.xs:sab.xe,sab.ys:sab.ye);
sab.nvol   =sab.vol.*sab.scope3d;
sab.ttvol  =sum(sum(sum(sab.nvol)));
sab.numgrd =length(find(sab.scope3d > 0.5));
sab.meanvol=sab.ttvol./sab.numgrd;
sab.grdct  =sab.meanvol./sab.vol;

%--- MAB ---%
mab.xs = 232; mab.xe = 402; 
mab.ys = 316; mab.ye = 465;
mab.vol    =vol(    :,mab.xs:mab.xe,mab.ys:mab.ye);
mab.funal  =funal(  :,mab.xs:mab.xe,mab.ys:mab.ye);
mab.scope3d=scope3d(:,mab.xs:mab.xe,mab.ys:mab.ye);
mab.nvol   =mab.vol.*mab.scope3d;
mab.ttvol  =sum(sum(sum(mab.nvol)));
mab.numgrd =length(find(mab.scope3d > 0.5));
mab.meanvol=mab.ttvol./mab.numgrd;
mab.grdct  =mab.meanvol./mab.vol;

%--- GOM ---%
gom.xs = 12;  gom.xe = 195;
gom.ys = 113; gom.ye = 257;
gom.vol    =vol(    :,gom.xs:gom.xe,gom.ys:gom.ye);
gom.funal  =funal(  :,gom.xs:gom.xe,gom.ys:gom.ye);
gom.funal(:,171:184,102:end)=0;
gom.scope3d=scope3d(:,gom.xs:gom.xe,gom.ys:gom.ye);
gom.scope3d(:,171:184,102:end)=0;
gom.nvol   =gom.vol.*gom.scope3d;
gom.ttvol  =sum(sum(sum(gom.nvol)));
gom.numgrd =length(find(gom.scope3d > 0.5));
gom.meanvol=gom.ttvol./gom.numgrd;
gom.grdct  =gom.meanvol./gom.vol;

%--- NIS ---%
nis.xs = 33 ; nis.xe = 95 ;
nis.ys = 146; nis.ye = 185;
nis.vol    =vol(    :,nis.ys:nis.ye,nis.xs:nis.xe);
nis.funal  =funal(  :,nis.ys:nis.ye,nis.xs:nis.xe);
nis.scope3d=scope3d(:,nis.ys:nis.ye,nis.xs:nis.xe);
nis.nvol   =nis.vol.*nis.scope3d;
nis.ttvol  =sum(sum(sum(nis.nvol)));
nis.numgrd =length(find(nis.scope3d > 0.5));
nis.meanvol=nis.ttvol./nis.numgrd;
nis.grdct  =nis.meanvol./nis.vol;

%=============================================================================%
%=                     Calculate various grid variables                      =%
%=============================================================================%
pdir='/Volumes/Black_box/Data/USeast-rtime/output/';
for i=6:10
  fname=[pdir,'useast_rtime_',sprintf('%04i',i),'.nc'];
  ncid=netcdf.open(fname,'nowrite');
    var=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'age_01'),...
                      [0 0 35 0],[402 482 1 720]);
  netcdf.close(ncid);
  if i==6
    rtime=squeeze(double(nanmean(var,4)))./5;
  else
    rtime=rtime+(squeeze(double(nanmean(var,4)))./5);
  end
end

sab.rtime=      rtime(    sab.xs:sab.xe,sab.ys:sab.ye).*sab.meanvol./...
          squeeze(vol(end,sab.xs:sab.xe,sab.ys:sab.ye));
mab.rtime=      rtime(    mab.xs:mab.xe,mab.ys:mab.ye).*mab.meanvol./...
          squeeze(vol(end,mab.xs:mab.xe,mab.ys:mab.ye));
gom.rtime=      rtime(    gom.xs:gom.xe,gom.ys:gom.ye).*gom.meanvol./...
          squeeze(vol(end,gom.xs:gom.xe,gom.ys:gom.ye));
nis.rtime=      rtime(    nis.xs:nis.xe,nis.ys:nis.ye).*nis.meanvol./...
          squeeze(vol(end,nis.xs:nis.xe,nis.ys:nis.ye));



