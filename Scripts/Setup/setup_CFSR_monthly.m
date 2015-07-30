%========================= setup_CFSR_monthly.m =========================%
%=                                                                      =%
%=  Written by Austin C. Todd, NCSU (2014) for personal use             =%
%=                                                                      =%
%=  Program reads in monthly averaged surface flux values from the CFSR =%
%=  Dataset and outputs them to new files to be used as forcing for the =%
%=  US East Coast climatology water age run. All CFSR variables are     =%
%=  interpolated to the US East coast ROMS grid for consistency (some   =%
%=  of the variables are available on ~0.3 degree grids, while others   =%
%=  are only available at 0.5 degree resolution). In addition to the    =%
%=  surface fluxes, the dQ/dSST field is also calculated and output to  =%
%=  file.                                                               =%
%=                                                                      =%
%========================================================================%

format long g; format compact;
addpath(genpath('/home/actodd/MYMATLAB'));

%========================================================================%
%=                  Set file names and directory paths                  =%
%========================================================================%
datadir  ='/raid0/Data/CFSR/Monthly/';

gridfile ='/raid0/actodd/USeast-age/Data/USeast-grid.nc';
fluxfile ='/raid0/actodd/USeast-age/Data/atmos/CFSR-fluxes.nc';
dsstfile ='/raid0/actodd/USeast-age/Data/atmos/CFSR-dQdSST.nc';
hdsstfile='/raid0/actodd/USeast-age/Data/atmos/CFSR-dQdSST-hycom.nc';
hycomfile='/raid0/Data/HYCOM/USeast-clim-monthly.nc';

%========================================================================%
%=                    Load in ROMS Grid information                     =%
%========================================================================%
ncid=netcdf.open(gridfile,'nowrite');
  roms.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho' )));
  roms.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho' )));
  roms.mask=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho')));
netcdf.close(ncid);
ncid=netcdf.open(hycomfile,'nowrite');
  roms.hycomsst=squeeze(double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,...
                               'temp'),[0 0 35 0],[402 482 1 12])));
  roms.hycomsst=permute(roms.hycomsst+273.15,[3 1 2]);
netcdf.close(ncid);

%========================================================================%
%=                   Load in Surface Momentum Fluxes                    =%
%========================================================================%
fname=[datadir,'stress/flxf06.gdas.U_FLX.SFC.grb2.nc'];

disp(['Reading U-stress from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  u.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  u.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  u.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  u.flux=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'U_FLX_L1_FcstAvg')));
netcdf.close(ncid);

fname=[datadir,'stress/flxf06.gdas.V_FLX.SFC.grb2.nc'];

disp(['Reading V-stress from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  v.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  v.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  v.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  v.flux=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'V_FLX_L1_FcstAvg')));
netcdf.close(ncid);

%========================================================================%
%=                  Load in winds and calculate speeds                  =%
%========================================================================%
fname=[datadir,'windspd/flxf06.gdas.WND.10m.grb2.nc'];

disp(['Reading Wind Speed from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  wnd.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  wnd.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  wnd.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  wnd.u   =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'U_GRD_L103_Avg')));
  wnd.v   =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'V_GRD_L103_Avg')));
netcdf.close(ncid);

wnd.flux=sqrt(wnd.u.^2 + wnd.v.^2);

%========================================================================%
%=                     Load in Downward Heat Flux                       =%
%========================================================================%
fname=[datadir,'heatflx/ocnh01.gdas.THFLX.SFC.grb2.nc'];

disp(['Reading Heat Flux from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  sh.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  sh.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  sh.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  sh.flux=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'THFLX_L1_FcstAvg')));
netcdf.close(ncid);

%========================================================================%
%=                     Load in Sensible Heat Flux                       =%
%========================================================================%
fname=[datadir,'sensible/flxf06.gdas.SHTFL.SFC.grb2.nc'];

disp(['Reading Sensible Heat Flux from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  snsbl.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  snsbl.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  snsbl.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  snsbl.flux=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'SHTFL_L1_Avg')));
netcdf.close(ncid);

%========================================================================%
%=                      Load in Latent Heat Flux                        =%
%========================================================================%
fname=[datadir,'latent/flxf06.gdas.LHTFL.SFC.grb2.nc'];

disp(['Reading Latent Heat Flux from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  ltnt.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  ltnt.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  ltnt.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  ltnt.flux=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'LHTFL_L1_Avg')));
netcdf.close(ncid);

%========================================================================%
%=                     Load in Evaporation-Precip                       =%
%========================================================================%
fname=[datadir,'eminusp/ocnh01.gdas.EMNP.SFC.grb2.nc'];

disp(['Reading E-P from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  emnp.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  emnp.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  emnp.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  emnp.flux=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'EMNP_L1_FcstAvg')));
netcdf.close(ncid);

%========================================================================%
%=              Load in Shortwave Radiation, calculate net              =%
%========================================================================%
%fname=[datadir,'pgbh01.gdas.DSWRF.SFC.grb2.nc'];
fname=[datadir,'rad/flxf06.gdas.DSWRF.SFC.grb2.nc'];

disp(['Reading Shortwave Radiation from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  swrad.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  swrad.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  swrad.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  swrad.down=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'DSWRF_L1_FcstAvg')));
netcdf.close(ncid);

fname=[datadir,'rad/flxf06.gdas.USWRF.SFC.grb2.nc'];
ncid=netcdf.open(fname,'nowrite');
  swrad.up=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'USWRF_L1_FcstAvg')));
netcdf.close(ncid);

%--- Make Net Shortwave Radiation (Down-Up) ---%
swrad.flux=swrad.down-swrad.up;

%========================================================================%
%=              Load in Longwave Radiation, calculate net               =%
%========================================================================%
fname=[datadir,'rad/flxf06.gdas.DLWRF.SFC.grb2.nc'];

disp(['Reading Shortwave Radiation from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  lwrad.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  lwrad.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  lwrad.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  lwrad.down=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'DLWRF_L1_Avg')));
netcdf.close(ncid);

fname=[datadir,'rad/flxf06.gdas.ULWRF.SFC.grb2.nc'];
ncid=netcdf.open(fname,'nowrite');
  lwrad.up=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ULWRF_L1_Avg')));
netcdf.close(ncid);

%--- Make Net Shortwave Radiation (Down-Up) ---%
lwrad.flux=lwrad.down-lwrad.up;

%========================================================================%
%=                      Load in 2m Air temperature                      =%
%========================================================================%
fname=[datadir,'temp/flxf06.gdas.TMP.2m.grb2.nc'];

disp(['Reading 2m Air temp from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  tmp2m.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  tmp2m.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  tmp2m.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  tmp2m.flux=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'TMP_L103_Avg')));
netcdf.close(ncid);

%========================================================================%
%=                      Load in Surface Temperature                     =%
%========================================================================%
fname=[datadir,'sst/flxf06.gdas.TMP.SFC.grb2.nc'];

disp(['Reading SST from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  sst.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  sst.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  sst.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  sst.flux=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'TMP_L1_Avg')));
netcdf.close(ncid);

%========================================================================%
%=                        Load in Surface Pressure                      =%
%========================================================================%
fname=[datadir,'pres/pgbh06.gdas.PRMSL.MSL.grb2.nc'];

disp(['Reading Pressure from ',fname]);
ncid=netcdf.open(fname,'nowrite');
  pres.time=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
  pres.lon =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon')));
  pres.lat =double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat')));
  pres.flux=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'PRMSL_L101_Avg')));
netcdf.close(ncid);

%========================================================================%
%=                                                                      =%
%= Loop through variables to calculate monthly means over entire record =%
%= then interpolate to the ROMS grid after masking out land points.     =% 
%=                                                                      =%
%========================================================================%
vars={'u';'v';'sh';'snsbl';'ltnt';'emnp';'swrad';'lwrad';'tmp2m';'sst';'pres';'wnd'};
    u.meanflux=zeros(size(    u.flux,1),size(    u.flux,2),12);
    v.meanflux=zeros(size(    v.flux,1),size(    v.flux,2),12);
   sh.meanflux=zeros(size(   sh.flux,1),size(   sh.flux,2),12);
snsbl.meanflux=zeros(size(snsbl.flux,1),size(snsbl.flux,2),12);
 ltnt.meanflux=zeros(size( ltnt.flux,1),size( ltnt.flux,2),12);
 emnp.meanflux=zeros(size( emnp.flux,1),size( emnp.flux,2),12);
swrad.meanflux=zeros(size(swrad.flux,1),size(swrad.flux,2),12);
lwrad.meanflux=zeros(size(lwrad.flux,1),size(lwrad.flux,2),12);

for var=1:length(vars)
  eval(['roms.',vars{var},'=zeros(12,size(roms.lon,1),size(roms.lon,2));']);

  eval([vars{var},'.flux(find(',vars{var},'.flux>1e35))=nan;']);
  eval(['temp=str2num(datestr(double(',vars{var},...
                                     '.time)./24+datenum(1979,1,1),5));']);
  for i=1:12
    a=find(temp==i);
    eval(['mnflx=squeeze(sum(',vars{var},'.flux(:,:,a),3))./length(a);']);
    eval(['[xlon,xlat]=meshgrid(',vars{var},'.lon-360,',vars{var},'.lat);']);

    %====================================================================%
    %=     Extrapolate the heat flux and E-P values to land points      =%
    %====================================================================%
    if (var==3 || var==4)
      flagGood=isfinite(mnflx);
      mnflx(~flagGood) = griddata(xlon( flagGood),...
                                  xlat( flagGood),...
      	                          mnflx(flagGood),...
                                  xlon(~flagGood),...
                                  xlat(~flagGood),'nearest');
    end   
    %====================================================================%
    %=                    Interpolate to ROMS grid                      =%
    %====================================================================%
    eval(['roms.',vars{var},'(i,:,:)=griddata(    xlon'',    xlat'',mnflx,',...
                                             'roms.lon,roms.lat);']);
    clear mnflx
  end
end
roms.time=[15.5,45,75.5,106,136.5,167,197.5,228.5,259,289.5,320,350.5];

%--- Adjust units on E-P ---%
roms.emnp=roms.emnp./100.*(1./24);
roms.htflx=roms.lwrad+roms.swrad-roms.snsbl-roms.ltnt;

%=========================================================================%
%=             Calculate saturated air specific humidity                 =%
%=========================================================================%
for i=1:size(roms.sst,1)
  for j=1:size(roms.sst,2)
    for k=1:size(roms.sst,3)
      roms.qs(     i,j,k)= (0.622./(0.01.*roms.pres(    i,j,k))).*...
                       10^(9.4051- (2553./roms.sst(     i,j,k)));
      roms.hycomqs(i,j,k)= (0.622./(0.01.*roms.pres(    i,j,k))).*...
                       10^(9.4051- (2553./roms.hycomsst(i,j,k)));
    end
  end
end

%=========================================================================%
%=                  Set Constants for calculations                       =%
%=========================================================================%
SBoltz=5.67e-8;     % Stefan-Boltzmann constant
Cp    =1.0048e3;   % Air specific heat at constant pressure
Ch    =1e-3;        % Bulk transfer coefficient for sensible heat
Ce    =1.15e-3;     % Bulk transfer coefficient for latent heat
L     =2.508e6;     % Latent heat of vaporization

%=========================================================================%
%=                       Calculate air density                           =%
%=========================================================================%
roms.rho= roms.pres./(roms.tmp2m.*287.058);   %rho = P/RT

%=========================================================================%
%=                         Calculate dQ/dSST                             =%
%=                                                                       =%
%=  Term 1: dQ/dSST due to Infrared radiation = 4*SBoltz*SST^3           =%
%=  Term 2: dQ/dSST due to Sensible heat flux = rho*Cp*Ch*Wind_Speed     =%
%=  Term 3: dQ/dSST due to Latent heat flux =                            =%
%=                           rho*Ce*L*Wind_Speed*2553*ln(10)*qs*(SST^-2) =%
%=                                                                       =%
%=========================================================================%
term1= 4.*SBoltz.*(roms.sst.^3);                        % term 1
term2= roms.rho.*(Cp*Ch).*roms.wnd;                         % term 2
term3= roms.rho.*(Ce*L*2553*log(10)).*roms.wnd.*roms.qs.*(roms.sst.^-2); % term 3

roms.dqdsst     =term1+term2+term3;
roms.hycomdqdsst=term1+term2+(roms.rho.*(Ce*L*2553*log(10)).*roms.wnd).*...
                             (roms.hycomqs.*(roms.hycomsst.^-2));

%========================================================================%
%=                     Output mean values to files                      =%
%========================================================================%
ncid=netcdf.open(fluxfile,'write');
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'time'),0,12,roms.time);
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'lon_rho'),roms.lon);
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'lat_rho'),roms.lat);
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'sustr' ),permute(roms.u,[2 3 1]));
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'svstr' ),permute(roms.v,[2 3 1]));
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'shflux'),permute(roms.htflx,[2 3 1]));
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'swflux'),permute(roms.emnp,[2 3 1]));
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'swrad' ),permute(roms.swrad,[2 3 1]));
netcdf.close(ncid);

ncid=netcdf.open(dsstfile,'write');
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'sst_time'),0,12,roms.time);
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'lon_rho'),roms.lon);
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'lat_rho'),roms.lat);
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'SST'),permute(roms.sst-273.15,[2 3 1]));
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'dQdSST'),permute(-roms.dqdsst,[2 3 1]));
netcdf.close(ncid);

ncid=netcdf.open(hdsstfile,'write');
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'sst_time'),0,12,roms.time);
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'lon_rho'),roms.lon);
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'lat_rho'),roms.lat);
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'SST'),permute(roms.hycomsst-273.15,[2 3 1]));
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'dQdSST'),permute(-roms.hycomdqdsst,[2 3 1]));
netcdf.close(ncid);






