%============================= plot_Obs.m ===============================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program the sampling locations from the NODC datasets.          =%
%=                                                                       =%
%=========================================================================%

%=========================================================================%
%=                    Set various paths and libraries                    =%
%=========================================================================%
addpath('/home/actodd/MYMATLAB/');
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));
addpath(genpath('/home/actodd/MYMATLAB/m_map/'));

%=========================================================================%
%=                Load ROMS netCDF file and read lat/lon                 =%
%=========================================================================%
roms.file=['/gpfs_share/actodd/USeast-age/output/2004/useast_his.nc'];
ncid=netcdf.open(roms.file,'nowrite');
  roms.lon        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho' ));
  roms.lat        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho' ));
  roms.h          =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'       ));
  roms.mask       =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'));
netcdf.close(ncid);

%=========================================================================%
%=                             Load CTD Data                             =%
%=========================================================================%
CTD.dir='/he_data/he/actodd/DATA/CTD/';
CTD.files=dir(CTD.dir);
ncid=netcdf.open([CTD.dir,CTD.files(3).name],'nowrite');
  CTD.lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon' ));
  CTD.lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat' ));
  CTD.time=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
  CTD.cast=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'cast'));
netcdf.close(ncid);

%=========================================================================%
%=                             Load MRB Data                             =%
%=========================================================================%
MRB.dir='/he_data/he/actodd/DATA/MRB/';
MRB.files=dir(MRB.dir);
ncid=netcdf.open([MRB.dir,MRB.files(3).name],'nowrite');
  MRB.lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon' ));
  MRB.lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat' ));
  MRB.time=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
  MRB.cast=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'cast'));
netcdf.close(ncid);

%=========================================================================%
%=                             Load XBT Data                             =%
%=========================================================================%
XBT.dir='/he_data/he/actodd/DATA/XBT/';
XBT.files=dir(XBT.dir);
ncid=netcdf.open([XBT.dir,XBT.files(3).name],'nowrite');
  XBT.lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon' ));
  XBT.lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat' ));
  XBT.time=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
  XBT.cast=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'cast'));
netcdf.close(ncid);

%=========================================================================%
%=                             Load PFL Data                             =%
%=========================================================================%
PFL.dir='/he_data/he/actodd/DATA/PFL/';
PFL.files=dir(PFL.dir);
ncid=netcdf.open([PFL.dir,PFL.files(3).name],'nowrite');
  PFL.lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon' ));
  PFL.lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat' ));
  PFL.time=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
  PFL.cast=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'cast'));
netcdf.close(ncid);

figure(1); clf;
m_proj('lambert','long',[min(min(roms.lon)) max(max(roms.lon))],...
                 'lat' ,[min(min(roms.lat)) max(max(roms.lat))]);
m_grid('tickdir','on');
hold on; 

m_contour(roms.lon,roms.lat,roms.mask,[0 0],'k');
m_text(-96.75,38,'NODC Observations','fontsize',14,'fontweight','bold');
m_text(-92.5,36.75,'2004-2006','fontsize',14,'fontweight','bold');
m_plot(XBT.lon,XBT.lat,'r.')
m_plot(CTD.lon,CTD.lat,'k.')
m_plot(MRB.lon,MRB.lat,'ms','markerfacecolor','m')
m_plot(PFL.lon,PFL.lat,'b.')
export_fig /home/actodd/PLOTS/USeast/validate/CTD_obs.png -png -r150 -painters
