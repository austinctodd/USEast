%======================== adjust_shelf_dQdSST.m ==========================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program is designed to read in calculated values of dQ/dSST     =%
%=  from the ROMS model setup, and to adjust values within shelf regions =%
%=  to be slightly higher than those calculated.  This is done in an     =%
%=  attempt to reduce excessive heating and cooling in these areas from  =%
%=  forcing with net heat fluxes.                                        =%
%=                                                                       =%
%=========================================================================%

%--- Add various libraries and paths ---%
addpath(genpath('/Users/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/Users/actodd/MYMATLAB/');
addpath(genpath('/Users/actodd/MYMATLAB/m_map'));
format long g; format compact;
set(0,'defaultfigurerenderer','zbuffer');

%=========================================================================%
%=                  Set file input/output directories                    =%
%=========================================================================%
plotdir='/Users/actodd/PLOTS/USeast/temp/frames/';
datadir='/Volumes/Black_box/Data/USeast/Data/atmos/';
climname='/Volumes/Black_box/Data/USeast/Data/clim/USeast-clim-monthly.nc';

%=========================================================================%
%=                      Load ROMS grid information                       =%
%=========================================================================%
disp('Reading in ROMS grid information');

fname=['/Volumes/Black_box/Data/USeast-age/output/clim/testdQdSST/',...
       'ERA/useast_his.nc'];

ncid=netcdf.open(fname,'nowrite');
 lon        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'    ));
 lat        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'    ));
 mask       =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'   ));
 h          =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'          ));
 shelf      =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_shelf'));
netcdf.close(ncid);

[sizelt,sizeln]=size(mask);
nanmask=mask./mask;

%=========================================================================%
%=                     Adjust shelf mask to gradient                     =%
%=========================================================================%
mask200 =mask'; mask200( find(h'>200 ))=0;

%--- Add some shelf waters in Western Gulf of Mexico ---%
mask200(118:119,39:41)=1; 
mask200(132:133,27)=1; 
mask200(131,28)=1;
mask200(151,19)=1;

%--- Fill in the keys ---%
mask200(185:188,174:188)=1; 
mask200(189:191,189:191)=1; 

%-- Extend shelf waters in W GOM ---%
mask200(118,42:43)=1;mask200(123,33)=1; mask200(123:124,32)=1;
mask200(127:128,30)=1;   mask200(129:131,29)=1;
mask200(132,28)=1; mask200(133,28)=1;
mask200(134,27)=1; mask200(135,26)=1;
mask200(148:150,20)=1; mask200(152,19)=1;
mask200(153,18)=1; 
mask200(120,37:39)=1; mask200(122,34)=1;
mask200(119,42)=1; mask200(118,44)=1;
mask200(151,20)=1;

%--- Gulf of Maine ---%
mask200(408:432,300:337)=1;
mask200(433:437,321:326)=1;
mask200(409:414,338:342)=1;
mask200(431:440,362:374)=1;
mask200(410:411,343)=1;
mask200(411,344)=1;

%--- Bermuda ---%
mask200(278:279,352:353)=0;

%get rid of central/south america mask
mask200(1:182,42:402)=0;

mask200=mask200';

for i=1:size(mask200,1)
  for j=1:size(mask200,2)
    if (mask(i,j)==0 && mask200(i,j)==0)
        mask200(i,j)=1;
    elseif (mask(i,j)==1 & mask200(i,j)==1)
        mask200(i,j)=mask200(i,j);
    elseif (mask(i,j)==1 && mask200(i,j)==0)
        mask200(i,j)=mask200(i,j);
    elseif (mask(i,j)==0 && mask200(i,j)==1)
        mask200(i,j)=mask200(i,j);
    end
  end
end

mask200=mask200';

%--- Carribean and S America ---%
mask200(52:220,197:402)=0;
mask200(118:180,141:199)=0;
mask200(89:94,124:135)=0;
mask200(41:50,307:402)=0;
mask200(277:281,351:355)=0;
mask200=mask200';

shelf=mask200;

for t=1:10
  for i=1:size(shelf,1)
      for j=1:size(shelf,2)
          shelf(i,j)=nanmean(nanmean(...
                            shelf(max([i-1 1]):min([i+1 size(shelf,1)]),...
                                  max([j-1 1]):min([j+1 size(shelf,2)]))));
      end
   end
end

%=========================================================================%
%=                         Load ERA temps/dQdSST                         =%
%=========================================================================%
ncid=netcdf.open([datadir,'ERA-dQdSST.nc'],'nowrite');
  era.lat   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
  era.lon   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
  era.dQdSST_raw=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dQdSST'));
  era.sst_raw   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'SST'));
netcdf.close(ncid);
era.lat   =double(era.lat);
era.lon   =double(era.lon);
era.dQdSST_raw=double(squeeze(era.dQdSST_raw));
era.sst_raw   =double(squeeze(era.sst_raw));

[x,y]=meshgrid(era.lon,era.lat);

%--- Interpolate to ROMS grid ---%
for i=1:12
    era.dQdSST(:,:,i)=griddata(x,y,squeeze(era.dQdSST_raw(:,:,i))',...
                               lon,lat);
    era.sst(   :,:,i)=griddata(x,y,squeeze(era.sst_raw(:,:,i))',...
                               lon,lat);
end

%=========================================================================%
%=                        Load CFSR temps/dQdSST                         =%
%=========================================================================%
ncid=netcdf.open([datadir,'CFSR-dQdSST.nc'],'nowrite');
  cfsr.lat   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'));
  cfsr.lon   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'));
  cfsr.dQdSST=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dQdSST'));
  cfsr.sst=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'SST'));
netcdf.close(ncid);
cfsr.lat   =double(cfsr.lat);
cfsr.lon   =double(cfsr.lon);
cfsr.dQdSST=double(squeeze(cfsr.dQdSST));
cfsr.sst   =double(squeeze(cfsr.sst   ));

%=========================================================================%
%=                      Load ERA/HYCOM temps/dQdSST                      =%
%=========================================================================%
ncid=netcdf.open([datadir,'ERA-dQdSST-hycom.nc'],'nowrite');
  era_hycom.lat   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
  era_hycom.lon   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
  era_hycom.dQdSST=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dQdSST'));
  era_hycom.sst=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'SST'));
netcdf.close(ncid);
era_hycom.lat   =double(era_hycom.lat);
era_hycom.lon   =double(era_hycom.lon);
era_hycom.dQdSST=double(squeeze(era_hycom.dQdSST));
era_hycom.sst   =double(squeeze(era_hycom.sst   ));

%=========================================================================%
%=                     Load CFSR/HYCOM temps/dQdSST                      =%
%=========================================================================%
ncid=netcdf.open([datadir,'CFSR-dQdSST-hycom.nc'],'nowrite');
  cfsr_hycom.lat   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'));
  cfsr_hycom.lon   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'));
  cfsr_hycom.dQdSST=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dQdSST'));
  cfsr_hycom.sst=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'SST'));
netcdf.close(ncid);
cfsr_hycom.lat   =double(cfsr_hycom.lat);
cfsr_hycom.lon   =double(cfsr_hycom.lon);
cfsr_hycom.dQdSST=double(squeeze(cfsr_hycom.dQdSST));
cfsr_hycom.sst   =double(squeeze(cfsr_hycom.sst   ));

%=========================================================================%
%=                    Write out adjusted dQdSST fields                   =%
%=========================================================================%
shelf_shift=zeros(size(lon,1),size(lon,2),12);
for i=1:12
  shelf_shift(:,:,i)=shelf.*0.25.*lat;
end

%--- ERA ---%
ncid=netcdf.open([datadir,'ERA-dQdSST-shelf.nc'],'write');
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'dQdSST'),...
                era.dQdSST-shelf_shift);
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'SST'),era.sst);
netcdf.close(ncid);

%--- ERA w/ HYCOM ---%
ncid=netcdf.open([datadir,'ERA-dQdSST-hycom-shelf.nc'],'write');
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'dQdSST'),...
                era_hycom.dQdSST-shelf_shift);
netcdf.close(ncid);

%--- CFSR ---%
ncid=netcdf.open([datadir,'CFSR-dQdSST-shelf.nc'],'write');
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'dQdSST'),...
                cfsr.dQdSST-shelf_shift);
netcdf.close(ncid);

%--- CFSR w/ HYCOM ---%
ncid=netcdf.open([datadir,'CFSR-dQdSST-hycom-shelf.nc'],'write');
  netcdf.putVar(ncid,netcdf.inqVarID(ncid,'dQdSST'),...
                cfsr_hycom.dQdSST-shelf_shift);
netcdf.close(ncid);






