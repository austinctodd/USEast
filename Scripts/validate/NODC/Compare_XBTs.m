%=========================== Compare_XBTs.m ==============================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program is designed to read in results from a the US East Coast =%
%=  water age ROMS simulation and compare to a whole set of XBT casts    =%
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
XBT.dir='/he_data/he/actodd/DATA/XBT/';
XBT.files=dir(XBT.dir);

%=========================================================================%
%=                  Get casts and contents of casts                      =%
%=========================================================================%
ncid=netcdf.open([XBT.dir,XBT.files(3).name],'nowrite');
  XBT.lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon' ));
  XBT.lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat' ));
  XBT.time=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
  XBT.cast=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'cast'));
netcdf.close(ncid);

%--- Adjust the time from the reference time ---%
XBT.time=XBT.time+datenum(1770,1,1);

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
  
stind=min(find(XBT.time<datenum(2005,1,1)));
enind=max(find(XBT.time<datenum(2005,1,1)));

for i=stind:enind
  
  disp(['Cast ',sprintf('%09i',XBT.cast(i))]);
    
  %--- Find closest lat/lon pairs on ROMS grid ---%
  %=======================================================================%
  %=              Find closest lat/lon points on ROMS grid               =%
  %=======================================================================%
  xmax=min(find(roms.lon(  : ,1)>=XBT.lon(i)));
  xmin=max(find(roms.lon(  : ,1)<=XBT.lon(i)));
  ymin=max(find(roms.lat(xmax,:)<=XBT.lat(i)));
  ymax=min(find(roms.lat(xmax,:)>=XBT.lat(i)));

  %=======================================================================%
  %=          Refine point on ROMS grid to find closest point            =%
  %=======================================================================%
  count=1; rdist=[]; xs=[]; ys=[];
  for x=max([1,xmin-3]):min([size(roms.lon,1),xmax+3])
    for y=max([1,ymin-3]):min([size(roms.lon,2),ymax+3])
      if roms.mask(x,y)
        rdist(count)=sw_dist([XBT.lat(i) roms.lat(x,y)],...
                             [XBT.lon(i) roms.lon(x,y)]);
        xs(count)=x; ys(count)=y;
      else
        rdist(count)=1e36;
        xs(count)=1; ys(count)=1;
      end
      count=count+1;  
    end
  end
  if (min(rdist)==1e36)
    disp(['WARNING: NO UNMASKED VALUES CAST',sprintf('%8i',XBT.cast(i))]);
    xind=1; yind=1;
  else
    a=find(rdist==min(rdist));
    if length(a)>1; a=a(1); end
    xind=xs(a); yind=ys(a);
  end
  clear dist xs ys x y a ymin ymax xmin xmax
  
  %--- Find closest ROMS time point ---%
  tind=max(find(roms.time<=XBT.time(i)));
  if (abs(XBT.time(i)-roms.time(tind+1))<abs(XBT.time(i)-roms.time(tind)))
      tind=tind+1;
  end 
  roms.casttime(i)=roms.time(tind);
  
  %=======================================================================%
  %=                           Get Cast data                             =%
  %=======================================================================%
  fname =[XBT.dir,'wod_',sprintf('%09i',XBT.cast(i)),'O.nc'];
  vnames=nc_vnames(fname);
  ncid=netcdf.open(fname,'nowrite');
  for j=1:length(vnames.Variables)
    if (strcmp(vnames.Variables(j).Name,'Temperature')==1)
      XBT.temp{  i}=double(netcdf.getVar(ncid,j-1));
    elseif (strcmp(vnames.Variables(j).Name,'Temperatuer_WODflag')==1)
      XBT.tflag{ i}=double(netcdf.getVar(ncid,j-1));
    elseif (strcmp(vnames.Variables(j).Name,'z')==1)
      XBT.z{     i}=double(netcdf.getVar(ncid,j-1));
    elseif (strcmp(vnames.Variables(j).Name,'WOD_cruise_identifier')==1)
      XBT.cruise{i}=     netcdf.getVar(ncid,j-1)';
    end
  end
  netcdf.close(ncid);
  
  %--- Do some Quality Control of the loaded dataset ---%
  if (length(XBT.temp)<i)
    XBT.temp{i}=[];
  else
    a=find(XBT.tflag{i}==0);
    if length(a>0)
      XBT.temp{i}=XBT.temp{i}(a);
      XBT.z{   i}=XBT.z{   i}(a);
    end
    
    %--- See if extraneous values exist ---%
    a=find(XBT.cruisetemp{i}{j}>-40 & XBT.temp{i}{j}<100);
    if (length(a)>0); 
      XBT.cruisetemp{i}{j}=XBT.cruisetemp{i}{j}(a);
      XBT.cruisez{   i}{j}=XBT.cruisez{   i}{j}(a);
    else
      XBT.temp{i}=[];
      XBT.z{   i}=[]
    end
  end
  
  %=======================================================================%
  %=                           Get ROMS data                             =%
  %=======================================================================%
  roms.casttime(i)=roms.time(tind);
  roms.castlon( i)=roms.lon(xind,yind);
  roms.castlat( i)=roms.lat(xind,yind);
  
  if (roms.mask(xind,yind)==0)
    roms.temp(i,:)=nan(1,36); roms.salt(i,:)=nan(1,36); 
    roms.z(   i,:)=nan(1,36);
    roms.tempi{i}=CTD.tz{i}.*nan; roms.salti{i}=CTD.sz{i}.*nan;
  else
    zeta=double(squeeze(netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'zeta'),...
                                    [xind-1 yind-1   tind-1],[1 1    1])));
    roms.temp(i,:)=double(squeeze(netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'temp'),...
                                    [xind-1 yind-1 0 tind-1],[1 1 36 1])));
    roms.z(i,:)=squeeze(set_depth(roms.Vtransform,roms.Vstretching,...
                                  roms.theta_s,roms.theta_b,roms.hc,...
                                  36,1,roms.h(xind,yind),zeta,0));
 
    %=====================================================================%
    %=                      Make interpolated casts                      =%
    %=====================================================================%
    if length(XBT.temp{i}>0)
      %--- Interpolate ROMS temp to Observed depths --%
      roms.tempi{i}=interp1(-roms.z(i,:),roms.temp(i,:),XBT.z{i});                   
    
      % Put into large array
      a=find(isnan(roms.tempi{i})==0);
      if (length(a)>0)
        if (i==stind)
          roms.temptot=roms.tempi{i}(a);
           XBT.temptot= XBT.temp{i}(a);
        else
          roms.temptot(end+1:end+length(a))=roms.tempi{i}(a);
           XBT.temptot(end+1:end+length(a))= XBT.temp{i}(a);
        end 
      end
    end
  end
end

%=========================================================================%
%=                         Group data by cruise                          =%
%=========================================================================%

%-- Find unique cruise names ---%
clear temp
for i=1:length(XBT.cruise)
  if length(XBT.cruise{i})>0
    temp(i,:)=XBT.cruise{i};
  end
end
cruisenames=unique(temp,'rows'); clear temp

%--- Initialize arrays for cruise casts ---%
roms.cruisetemp=[]; roms.cruisetempi=[]; 
 XBT.cruisetemp=[];  
 
%--- All data along one array (for cruise stats) ---%
roms.cruisetemptot   =[]; roms.cruisetemptot{1}=[]; 
 XBT.cruisetemptot   =[];  XBT.cruisetemptot{1}=[];  
 
for i=1:size(cruisenames,1)
  a=find(strcmp(XBT.cruise,cruisenames(i,:))==1);
  
  %--- Loop through each cast from a cruise ---%
  for j=1:length(a)
      
    %--- Sort ROMS data into cruise and cast ---%
    roms.cruisetemp{   i}{j}=roms.temp( a(j),:);
    roms.cruisetempi{  i}{j}=roms.tempi{a(j)  };
    roms.cruisetime(   i, j)=roms.casttime(a(j));
    roms.cruiselon(    i, j)=roms.castlon( a(j));
    roms.cruiselat(    i, j)=roms.castlat( a(j));
    roms.cruisez{      i}{j}=roms.z(       a(j),:);
    
    %--- Now for CTD data ---%
    XBT.cruisetemp{   i}{j}= XBT.temp{ a(j)};
    XBT.cruisetime(   i, j)= XBT.time( a(j));
    XBT.cruiselon(    i, j)= XBT.lon(  a(j));
    XBT.cruiselat(    i, j)= XBT.lat(  a(j));   
    XBT.cruisez{      i}{j}= XBT.z(    a(j));   

    %--- Make long array for entire cruise ---%
    b=find(isnan(roms.tempi{a(j)})==0);
    if (length(roms.cruisetemptot)<i && length(b)>0)
      roms.cruisetemptot{i}   =roms.tempi{a(j)}(b);
       XBT.cruisetemptot{i}= XBT.temp{a(j)}(b);
    end
    if (length(roms.cruisetemptot)==i && length(b)>0)
      roms.cruisetemptot{i}(end+1:end+length(b))=roms.tempi{a(j)}(b);
       XBT.cruisetemptot{i}(end+1:end+length(b))= XBT.temp{a(j)}(b);
    end
  end
  %=======================================================================%
  %               Now calculate some stats for each cruise               =%
  %=======================================================================%
 
  %--- RMS ---%
  Tcruise.rms(i)=sqrt(mean((XBT.cruisetemptot{i}-...
                           roms.cruisetemptot{i}).^2));
  %--- R/P ---%
  [r,p]= corrcoef(XBT.cruisetemptot{i},roms.cruisetemptot{i});
  Tcruise.r(i)=r(1,2); Tcruise.p(i)=p(1,2);
  
  %--- Standard Deviation ---%
  Tcruise.std(1,i)=std(roms.cruisetemptot{i});
  Tcruise.std(2,i)=std( XBT.cruisetemptot{i});
  
end

netcdf.close(ncid4); netcdf.close(ncid5);

%save XBT_ROMS_data.mat XBT roms -MAT
