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
addpath(genpath('/home/actodd/MYMATLAB_temp/ROMS-matlab/'));
%addpath('/home/actodd/MYMATLAB/');
%addpath(genpath('/home/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=========================================================================%
%=                  Set file input/output directories                    =%
%=========================================================================%
plotdir='/he_data/he/actodd/PLOTS/USeast/validate/KE/';
datadir='/gpfs_common/he_share/actodd/USeast-age/output/clim/fwdrun/';
matdir ='/he_data/he/actodd/DATA/eke/';

%=========================================================================%
%=                      Load ROMS grid information                       =%
%=========================================================================%
disp('Reading in ROMS output');

ncid=netcdf.open([datadir,'part1/useast_his.nc'],'nowrite');
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
templat=zeros(2922*6,201)

ncid1=netcdf.open([datadir,'part1/useast_his.nc'],'nowrite');
ncid2=netcdf.open([datadir,'part2/useast_his.nc'],'nowrite');

for tdim=0:2921

  disp(['Day ',sprintf('%7.3f',(tdim/8)+1)]);

  for yr=1:6

    %--- Load temp and SSH, calculate depth ---%
    inds=[0+tdim,2922+tdim,0+tdim,2922+tdim,5844+tdim,8766+tdim];
    if (yr <3) ncid=ncid1; else ncid=ncid2; end
    temp=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'temp'),...
                         [201 289 0 inds(yr)],[201 111 36 1]);
    zeta=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'),...
                         [201 289   inds(yr)],[201 111    1]);
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

          templat(count(i),i)=lat2(i,x(1));
          tempval(count(i),i)=temp2(x(1));
	  count(i)=count(i)+1;
        end
      end
    end
  end
end
netcdf.close(ncid);

save /gpfs_share/actodd/GS_positions.mat templat tempval count -MAT

