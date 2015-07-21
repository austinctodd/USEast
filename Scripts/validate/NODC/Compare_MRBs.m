%=========================== Compare_MRBs.m ==============================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program is designed to read in results from a the US East Coast =%
%=  water age ROMS simulation and compare to a whole set of moored buoys =%
%=  from the NODC datasets.                                              =%
%=                                                                       =%
%=========================================================================%

%=========================================================================%
%=                    Set various paths and libraries                    =%
%=========================================================================%
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));

%=========================================================================%
%=           Define data directory and get contents of directory         =%
%=========================================================================%
MRB.dir='/he_data/he/actodd/DATA/MRB/';
MRB.files=dir(MRB.dir);

%=========================================================================%
%=                  Get casts and contents of casts                      =%
%=========================================================================%
ncid=netcdf.open([MRB.dir,MRB.files(3).name],'nowrite');
  MRB.lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon' ));
  MRB.lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat' ));
  MRB.time=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
  MRB.cast=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'cast'));
netcdf.close(ncid);

%--- Adjust the time from the reference time ---%
MRB.time=MRB.time+datenum(1770,1,1);

%=========================================================================%
%=               Load ROMS netCDF files and read lat/lon                 =%
%=========================================================================%
roms.file2004=['/gpfs_share/actodd/USeast-age/output/2004/useast_his.nc'];
roms.file2005=['/gpfs_share/rhe/actodd/USeast-age/output/OBC/2005/',...
               'useast_his.nc'];
ncid4=netcdf.open(roms.file2004,'nowrite');
ncid5=netcdf.open(roms.file2005,'nowrite');

roms.lon        =netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'lon_rho' ));
roms.lat        =netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'lat_rho' ));
roms.h          =netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'h'       ));
roms.mask       =netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'mask_rho'));
roms.Vtransform =netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'Vtransform' ));
roms.Vstretching=netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'Vstretching'));
roms.theta_s    =netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'theta_s'    ));
roms.theta_b    =netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'theta_b'    ));
roms.hc         =netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'hc'         ));
         
%=========================================================================%
%=                            Get cast data                              =%
%=========================================================================%

%---- First, go through 2004 data ---%
roms.time=netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'ocean_time'));
roms.time=roms.time./86400+datenum(2004,1,1);
  
stind=min(find(MRB.time<datenum(2005,1,1)));
enind=max(find(MRB.time<datenum(2005,1,1)));

for i=stind:enind
  
  disp(['Cast ',sprintf('%09i',MRB.cast(i))]);
    
  %--- Find closest lat/lon pairs on ROMS grid ---%
  %=======================================================================%
  %=              Find closest lat/lon points on ROMS grid               =%
  %=======================================================================%
  xmax=min(find(roms.lon(  : ,1)>=MRB.lon(i)));
  xmin=max(find(roms.lon(  : ,1)<=MRB.lon(i)));
  ymin=max(find(roms.lat(xmax,:)<=MRB.lat(i)));
  ymax=min(find(roms.lat(xmax,:)>=MRB.lat(i)));

  %=======================================================================%
  %=          Refine point on ROMS grid to find closest point            =%
  %=======================================================================%
  count=1; rdist=[]; xs=[]; ys=[];
  for x=max([1,xmin-3]):min([size(roms.lon,1),xmax+3])
    for y=max([1,ymin-3]):min([size(roms.lon,2),ymax+3])
      if roms.mask(x,y)
        rdist(count)=sw_dist([MRB.lat(i) roms.lat(x,y)],...
                             [MRB.lon(i) roms.lon(x,y)]);
        xs(count)=x; ys(count)=y;
      else
        rdist(count)=1e36;
        xs(count)=1; ys(count)=1;
      end
      count=count+1;  
    end
  end
  if (min(rdist)==1e36)
    disp(['WARNING: NO UNMASKED VALUES CAST',sprintf('%8i',MRB.cast(i))]);
    xind=1; yind=1;
  else
    a=find(rdist==min(rdist));
    xind=xs(a); yind=ys(a);
  end
  clear dist xs ys x y a ymin ymax xmin xmax
  
  %--- Find closest ROMS time point ---%
  tind=max(find(roms.time<=MRB.time(i)));
  if (abs(MRB.time(i)-roms.time(tind+1))<abs(MRB.time(i)-roms.time(tind)))
      tind=tind+1;
  end 
  
  %=======================================================================%
  %=                           Get Cast data                             =%
  %=======================================================================%
  ncid=netcdf.open([MRB.dir,'wod_',sprintf('%09i',MRB.cast(i)),'O.nc'],...
                   'nowrite');
    MRB.temp{i}=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Temperature')));
    %CTD.salt{i}=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Salinity'   )));
    MRB.z{   i}=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'z'          )));
  netcdf.close(ncid);
  
  %=======================================================================%
  %=                           Get ROMS data                             =%
  %=======================================================================%
  
  zeta=double(squeeze(netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'zeta'),...
                                    [xind-1 yind-1   tind-1],[1 1    1])));
  temp=double(squeeze(netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'temp'),...
                                    [xind-1 yind-1 0 tind-1],[1 1 36 1])));
 % salt=double(squeeze(netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'salt'),...
  %                                  [xind-1 yind-1 0 tind-1],[1 1 36 1])));

  z=squeeze(set_depth(roms.Vtransform,roms.Vstretching,...
                       roms.theta_s,roms.theta_b,roms.hc,...
                       36,1,roms.h(xind,yind),zeta,0));
  
  roms.temp{i}=interp1(-z,temp,MRB.z{i});                   
  %roms.salt{i}=interp1(-z,salt,CTD.z{i});                   
                  
  
end










netcdf.close(ncid4); netcdf.close(ncid5);