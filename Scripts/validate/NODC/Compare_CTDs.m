%=========================== Compare_CTDs.m ==============================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program is designed to read in results from a the US East Coast =%
%=  water age ROMS simulation and compare to a whole set of CTD casts    =%
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
CTD.dir='/he_data/he/actodd/DATA/CTD/';
CTD.files=dir(CTD.dir);

%=========================================================================%
%=                  Get casts and contents of casts                      =%
%=========================================================================%
ncid=netcdf.open([CTD.dir,CTD.files(3).name],'nowrite');
  CTD.lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon' ));
  CTD.lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat' ));
  CTD.time=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
  CTD.cast=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'cast'));
netcdf.close(ncid);

%--- Adjust the time from the reference time ---%
CTD.time=CTD.time+datenum(1770,1,1);

%=========================================================================%
%=               Load ROMS netCDF files and read lat/lon                 =%
%=========================================================================%
roms.file2004=['/gpfs_share/actodd/USeast-age/output/2004/useast_his.nc'];
roms.file2005=['/gpfs_share/actodd/USeast-age/output/2005/useast_his.nc'];
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
roms.time4=netcdf.getVar(ncid4,netcdf.inqVarID(ncid4,'ocean_time'));
roms.time4=roms.time4./86400+datenum(2004,1,1);

%---- Add in 2005-2006 data ---%
roms.time5=netcdf.getVar(ncid5,netcdf.inqVarID(ncid5,'ocean_time'));
roms.time5=roms.time5./86400+datenum(2004,1,1);
clear tmp;

stind=min(find(CTD.time>=datenum(2004,1,1)));
enind=max(find(CTD.time<datenum(2007,1,1)));

for i=stind:enind
  
  disp(['Cast ',sprintf('%09i',CTD.cast(i))]);
    
  %--- Find closest lat/lon pairs on ROMS grid ---%
  %=======================================================================%
  %=              Find closest lat/lon points on ROMS grid               =%
  %=======================================================================%
  xmax=min(find(roms.lon(  : ,1)>=CTD.lon(i)));
  xmin=max(find(roms.lon(  : ,1)<=CTD.lon(i)));
  ymin=max(find(roms.lat(xmax,:)<=CTD.lat(i)));
  ymax=min(find(roms.lat(xmax,:)>=CTD.lat(i)));

  %=======================================================================%
  %=          Refine point on ROMS grid to find closest point            =%
  %=======================================================================%
  count=1; rdist=[]; xs=[]; ys=[];
  for x=max([1,xmin-3]):min([size(roms.lon,1),xmax+3])
    for y=max([1,ymin-3]):min([size(roms.lon,2),ymax+3])
      if roms.mask(x,y)
        rdist(count)=sw_dist([CTD.lat(i) roms.lat(x,y)],...
                             [CTD.lon(i) roms.lon(x,y)]);
        xs(count)=x; ys(count)=y;
      else
        rdist(count)=1e36;
        xs(count)=1; ys(count)=1;
      end
      count=count+1;  
    end
  end
  if (min(rdist)==1e36)
    disp(['WARNING: NO UNMASKED VALUES CAST',sprintf('%8i',CTD.cast(i))]);
    xind=1; yind=1;
  else
    a=find(rdist==min(rdist));
    if length(a)>1; a=a(1); end
    xind=xs(a); yind=ys(a);
  end
  clear dist xs ys x y a ymin ymax xmin xmax
  
  %--- Find closest ROMS time point ---%
  if (CTD.time(i)>=datenum(2005,1,1))
    tind=max(find(roms.time5<=CTD.time(i)));
    if (abs(CTD.time(i)-roms.time5(tind+1))<abs(CTD.time(i)-roms.time5(tind)))
        tind=tind+1;
    end 
    roms.casttime(i)=roms.time5(tind);
    ncind=ncid5;
  else
    tind=max(find(roms.time4<=CTD.time(i)));
    if (abs(CTD.time(i)-roms.time4(tind+1))<abs(CTD.time(i)-roms.time4(tind)))
        tind=tind+1;
    end 
    roms.casttime(i)=roms.time4(tind);
    ncind=ncid4;
  end
  
  %=======================================================================%
  %=                           Get Cast data                             =%
  %=======================================================================%
  fname =[CTD.dir,'wod_',sprintf('%09i',CTD.cast(i)),'O.nc'];
  vnames=nc_vnames(fname);
  ncid=netcdf.open(fname,'nowrite');
  for j=1:length(vnames.Variables)
    if (strcmp(vnames.Variables(j).Name,'Temperature')==1)
      CTD.temp{  i}=double(netcdf.getVar(ncid,j-1));
    elseif (strcmp(vnames.Variables(j).Name,'Temperature_WODflag')==1)
      CTD.tflag{ i}=double(netcdf.getVar(ncid,j-1));
    elseif (strcmp(vnames.Variables(j).Name,'Salinity')==1)
      CTD.salt{  i}=double(netcdf.getVar(ncid,j-1));
    elseif (strcmp(vnames.Variables(j).Name,'Salinity_WODflag')==1)
      CTD.sflag{ i}=double(netcdf.getVar(ncid,j-1));
    elseif (strcmp(vnames.Variables(j).Name,'z')==1)
      CTD.tz{     i}=double(netcdf.getVar(ncid,j-1));
      CTD.sz{     i}=double(netcdf.getVar(ncid,j-1));
    elseif (strcmp(vnames.Variables(j).Name,'WOD_cruise_identifier')==1)
      CTD.cruise{i}=     netcdf.getVar(ncid,j-1)';
    end
  end
  netcdf.close(ncid);
  
  %--- Do some Quality Control of the loaded dataset ---%
  if (length(CTD.temp)<i)
    CTD.temp{i}=[];
  else
    a=find(CTD.tflag{i}==0);
    if length(a>0)
      CTD.temp{i}=CTD.temp{i}(a);
      CTD.tz{  i}=CTD.tz{  i}(a);
    end
    %--- See if extraneous values exist ---%
    a=find(CTD.temp{i}>-40 & CTD.temp{i}<100);
    if (length(a)>0); 
      CTD.temp{i}=CTD.temp{i}(a);      
      CTD.tz{  i}=CTD.tz{  i}(a);
    else
      CTD.temp{i}=[];
      CTD.tz{  i}=[];
    end
  end
  if (length(CTD.salt)<i)
    CTD.salt{i}=[];
  else
    a=find(CTD.sflag{i}==0);
    if length(a>0)
      CTD.salt{i}=CTD.salt{i}(a);
      CTD.sz{  i}=CTD.sz{  i}(a);
    end
    %--- See if extraneous values exist ---%
    a=find(CTD.salt{i}>-40 & CTD.salt{i}<50);
    if (length(a)>0); 
      CTD.salt{i}=CTD.salt{i}(a);
      CTD.sz{  i}=CTD.sz{  i}(a);
    else
      CTD.salt{i}=[];
      CTD.sz{  i}=[];
    end
  end
  
  %=======================================================================%
  %=                           Get ROMS data                             =%
  %=======================================================================%
  if (CTD.time(i)<datenum(2005,1,1))
    roms.casttime(i)=roms.time4(tind);
  else
    roms.casttime(i)=roms.time5(tind);
  end
  roms.castlon( i)=roms.lon(xind,yind);
  roms.castlat( i)=roms.lat(xind,yind);
  
  if (roms.mask(xind,yind)==0)
    roms.temp(i,:)=nan(1,36); roms.salt(i,:)=nan(1,36); 
    roms.z(   i,:)=nan(1,36);
    roms.tempi{i}=CTD.tz{i}.*nan; roms.salti{i}=CTD.sz{i}.*nan;
  else
    zeta=double(squeeze(netcdf.getVar(ncind,netcdf.inqVarID(ncind,'zeta'),...
                                    [xind-1 yind-1   tind-1],[1 1    1])));
    roms.temp(i,:)=double(squeeze(netcdf.getVar(ncind,netcdf.inqVarID(ncind,'temp'),...
                                    [xind-1 yind-1 0 tind-1],[1 1 36 1])));
    roms.salt(i,:)=double(squeeze(netcdf.getVar(ncind,netcdf.inqVarID(ncind,'salt'),...
                                    [xind-1 yind-1 0 tind-1],[1 1 36 1])));

    roms.z(i,:)=squeeze(set_depth(roms.Vtransform,roms.Vstretching,...
                                  roms.theta_s,roms.theta_b,roms.hc,...
                                  36,1,roms.h(xind,yind),zeta,0));

  %=======================================================================%
  %=                       Make interpolated casts                       =%
  %=======================================================================%
    if length(CTD.temp{i}>0)
      %--- Interpolate ROMS temp to Observed depths --%
      roms.tempi{i}=interp1(-roms.z(i,:),roms.temp(i,:),CTD.tz{i});                   
    
      % Put into large array
      a=find(isnan(roms.tempi{i})==0);
      if (length(a)>0)
        if (i==stind)
          roms.temptot=roms.tempi{i}(a);
           CTD.temptot= CTD.temp{i}(a);
        else
          roms.temptot(end+1:end+length(a))=roms.tempi{i}(a);
           CTD.temptot(end+1:end+length(a))= CTD.temp{i}(a);
        end 
      end
    end
  
    if length(CTD.salt{i}>0)  
      %--- Interpolate ROMS salinity to Observed depths ---%
      roms.salti{i}=interp1(-roms.z(i,:),roms.salt(i,:),CTD.sz{i});                   

      % Put into large array
      a=find(isnan(roms.salti{i})==0);
      if (length(a)>0)
        if (i==stind)
          roms.salttot=roms.salti{i}(a);
           CTD.salttot= CTD.salt{ i}(a);
        else
          roms.salttot(end+1:end+length(a))=roms.salti{i}(a);
           CTD.salttot(end+1:end+length(a))= CTD.salt{ i}(a);
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
for i=1:length(CTD.cruise)
  if length(CTD.cruise{i})>0
    temp(i,:)=CTD.cruise{i};
  end
end
cruisenames=unique(temp,'rows'); clear temp

%--- Initialize arrays for cruise casts ---%
roms.cruisetemp=[]; roms.cruisetempi=[]; 
roms.cruisesalt=[]; roms.cruisesalti=[];
 CTD.cruisetemp=[];  CTD.cruisesalt =[];
 
%--- All data along one array (for cruise stats) ---%
roms.cruisetemptot   =[]; roms.cruisesalttot   =[];
roms.cruisetemptot{1}=[]; roms.cruisesalttot{1}=[];
 CTD.cruisetemptot   =[];  CTD.cruisesalttot   =[];
 CTD.cruisetemptot{1}=[];  CTD.cruisesalttot{1}=[];
 
for i=1:size(cruisenames,1)
  a=find(strcmp(CTD.cruise,cruisenames(i,:))==1);
  
  %--- Loop through each cast from a cruise ---%
  for j=1:length(a)
      
    %--- Sort ROMS data into cruise and cast ---%
    roms.cruisesalt{   i}{j}=roms.salt( a(j),:);
    roms.cruisesalti{  i}{j}=roms.salti{a(j)  };
    roms.cruisetemp{   i}{j}=roms.temp( a(j),:);
    roms.cruisetempi{  i}{j}=roms.tempi{a(j)  };
    roms.cruisetime(   i, j)=roms.casttime(a(j));
    roms.cruiselon(    i, j)=roms.castlon( a(j));
    roms.cruiselat(    i, j)=roms.castlat( a(j));
    roms.cruisez{      i}{j}=roms.z(       a(j),:);
    
    %--- Now for CTD data ---%
    CTD.cruisetemp{   i}{j}= CTD.temp{a(j)};
    CTD.cruisesalt{   i}{j}= CTD.salt{a(j)};
    CTD.cruisetime(   i, j)= CTD.time(a(j));
    CTD.cruiselon(    i, j)= CTD.lon( a(j));
    CTD.cruiselat(    i, j)= CTD.lat( a(j));   
    CTD.cruisetz{     i}{j}= CTD.tz(  a(j));   
    CTD.cruisesz{     i}{j}= CTD.sz(  a(j));   

    %--- Make long array for entire cruise ---%
    b=find(isnan(roms.tempi{a(j)})==0);
    c=find(isnan(roms.salti{a(j)})==0);
    if (length(roms.cruisetemptot)<i && length(b)>0)
      roms.cruisetemptot{i}   =roms.tempi{a(j)}(b);
       CTD.cruisetemptot{i}= CTD.temp{a(j)}(b);
    end
    if (length(roms.cruisesalttot)<i && length(c)>0)
      roms.cruisesalttot{i}=roms.salti{a(j)}(c);
       CTD.cruisesalttot{i}= CTD.salt{a(j)}(c);
    end      
    if (length(roms.cruisetemptot)==i && length(b)>0)
      roms.cruisetemptot{i}(end+1:end+length(b))=roms.tempi{a(j)}(b);
       CTD.cruisetemptot{i}(end+1:end+length(b))= CTD.temp{a(j)}(b);
    end
    if (length(roms.cruisesalttot)==i && length(c)>0)
      roms.cruisesalttot{i}(end+1:end+length(c))=roms.salti{a(j)}(c);
       CTD.cruisesalttot{i}(end+1:end+length(c))= CTD.salt{a(j)}(c);
    end      
  end
  %=======================================================================%
  %               Now calculate some stats for each cruise               =%
  %=======================================================================%
 
  if (length(roms.cruisetemptot)==i)
    %--- RMS ---%
    Tcruise.rms(i)=sqrt(nanmean((         CTD.cruisetemptot{i} -...
                                 nanmean( CTD.cruisetemptot{i})-...
                                         roms.cruisetemptot{i} -...
                                nanmean(roms.cruisetemptot{i})).^2));
    Scruise.rms(i)=sqrt(nanmean((         CTD.cruisesalttot{i} -...
                                 nanmean( CTD.cruisesalttot{i})-...
                                         roms.cruisesalttot{i} -...
                                 nanmean(roms.cruisesalttot{i})).^2));

    %--- Standard Deviation ---%
    Tcruise.std(1,i)=nanstd(roms.cruisetemptot{i});
    Tcruise.std(2,i)=nanstd( CTD.cruisetemptot{i});
    Scruise.std(1,i)=nanstd(roms.cruisesalttot{i});
    Scruise.std(2,i)=nanstd( CTD.cruisesalttot{i});
  
    %--- R/P ---%
    %[r,p]= corrcoef(CTD.cruisetemptot{i},roms.cruisetemptot{i});
    %Tcruise.r(i)=r(1,2); Tcruise.p(i)=p(1,2);    
    Tcruise.r(i)=nanmean((CTD.cruisetemptot{i}-nanmean( CTD.cruisetemptot{i}).*...
                        (roms.cruisetemptot{i}-nanmean(roms.cruisetemptot{i})))./...
                 (Tcruise.std(1,i)*Tcruise.std(2,i)));
    
    %[r,p]= corrcoef(CTD.cruisesalttot{i},roms.cruisesalttot{i});
    %Scruise.r(i)=r(1,2); Scruise.p(i)=p(1,2);
    Scruise.r(i)=nanmean((CTD.cruisesalttot{i}-nanmean( CTD.cruisesalttot{i}).*...
                        (roms.cruisesalttot{i}-nanmean(roms.cruisesalttot{i})))./...
                 (Scruise.std(1,i)*Scruise.std(2,i)));
  
  else
    Tcruise.rms(  i)=nan; Scruise.rms(  i)=nan;
    Tcruise.r(    i)=nan; Scruise.r(    i)=nan;
    Tcruise.std(:,i)=nan; Scruise.std(:,i)=nan;
  end
end

netcdf.close(ncid4); netcdf.close(ncid5);

save CTD_ROMS_data.mat CTD roms -MAT
