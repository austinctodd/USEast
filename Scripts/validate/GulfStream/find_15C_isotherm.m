%========================= find_15C_isotherm.m ===========================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This reads the temperature field from the US East water age model    =%
%=  to find the location of the 15 degree isotherm at 200m depth (which  =%
%=  is indicative of the Gulf Stream Wall).  The position is extracted   =%
%=  in the form of a latitude at each grid longitude (using interpolation=%
%=  routines.                                                            =%
%=                                                                       =%
%=========================================================================%

%--- Add various libraries and paths ---%
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/home/actodd/MYMATLAB/');
addpath(genpath('/home/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=========================================================================%
%=                  Set file input/output directories                    =%
%=========================================================================%
plotdir='/he_data/he/actodd/PLOTS/USeast/validate/KE/';
datadir='/gpfs_common/he_share/actodd/USeast-age/output/clim/';
matdir ='/he_data/he/actodd/DATA/eke/';

%=========================================================================%
%=                      Load ROMS grid information                       =%
%=========================================================================%
disp('Reading in ROMS output');

ncid=netcdf.open([datadir,'useast_his_0001.nc'],'nowrite');
 Vstretching=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vstretching'));
 Vtransform =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vtransform' ));
 theta_s    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_s'    ));
 theta_b    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_b'    ));
 hc         =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hc'         ));
 lon        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'    ));
 lat        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'    ));
 mask       =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'   ));
 h          =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'          ));
netcdf.close(ncid);

%--- Set up depth lat/lon points ---%
lon2=lon(202:402,290:400);
lat2=lat(202:402,290:400);
h2  =h(  202:402,290:400);

%=========================================================================%
%=                       Loop through time record                        =%
%=========================================================================%

count=ones(1,201);

fname=[datadir,'daily_avg.nc'];
ncid=netcdf.open(fname,'nowrite');
for j=1:364

  disp(['Day ',sprintf('%03i',j)]);

  %--- Load temp and SSH, calculate depth ---%
  temp=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'temp'),...
                       [201 289 0 j-1],[201 111 36 1]);
  zeta=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'),...
                       [201 289   j-1],[201 111    1]);
  z=set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,36,1,...
              h(202:402,290:400),zeta,0);

  %--- Go along longitude direction, find where 15C occurs ---% 
  for i=1:201
    for ii=1:111
      if(h2(i,ii)<200)
        temp2(ii)=nan;
      else
        temp2(ii)=interp1(squeeze(squeeze(   z(i,ii,:))),...
                          squeeze(squeeze(temp(i,ii,:))),-200);
      end
    end
    a=find(isnan(temp2)==0);
    if length(a)>0
      x=find( abs(temp2-15)==nanmin(abs(temp2-15)));

      if (count(i)==1)
        templat(j,1)=lat2(i,x(1));
	tempval{j,1)=temp2(x(1));
      else      
        templat(j,i)=lat2(i,x(1));
        tempval(j,i)=temp2(x(1));
	count(i)=count(i)+1;
      end
    end
  end
end
netcdf.close(ncid);

save /he_data/he/actodd/DATA/GS_position/avg_lats.mat templat tempval count -MAT

