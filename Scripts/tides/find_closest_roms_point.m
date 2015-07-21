function [xind,yind]=find_closest_roms_point(roms,ss,stnID);
%FIND_CLOSEST_ROMS_POINT
%=========================================================================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  Program finds the ROMS grid lat/lon pairs that are closest to a tide =%
%=  gauge station given by ID stnID                                      =%
%=                                                                       =%
%=  INPUT VARIABLES:                                                     =%
%=    roms  ::  structured array with variables roms.lon and roms.lat    =%
%=    ss    ::  structured array read in from US_TideInfo.mat            =%
%=    stnID ::  station identifier to find closest lat/lon pair for      =%
%=                                                                       =%
%=  OUTPUT VARIABLES:                                                    =%
%=    xind  ::  i-index for roms grid to extract data                    =%
%=    yind  ::  j-index for roms grid to extract data                    =%
%=                                                                       =%
%=========================================================================%

f_ind=find(ss.ID==stnID);
if isempty(f_ind)
    error(['cannot find station identifier ',stnID]);
end

%--- Find the closest lat/lon pairs for station ---%
xmax=min(find(roms.lon(1,   :)>=ss.lon(f_ind)));
xmin=max(find(roms.lon(1,   :)<=ss.lon(f_ind)));
ymin=max(find(roms.lat(:,xmax)<=ss.lat(f_ind)));
ymax=min(find(roms.lat(:,xmax)>=ss.lat(f_ind)));

count=1; rdist=[]; xs=[]; ys=[];
for x=max([1,xmin-3]):min([size(roms.lon,2),xmax+3])
  for y=max([1,ymin-3]):min([size(roms.lon,1),ymax+3])
    if roms.mask(y,x)
      rdist(count)=sw_dist([ss.lat(f_ind) roms.lat(y,x)],...
                           [ss.lon(f_ind) roms.lon(y,x)]);
      xs(count)=x; ys(count)=y;
    else
      rdist(count)=1e36;
      xs(count)=1; ys(count)=1;
    end
    count=count+1;  
  end
end
if (min(rdist)==1e36)
  disp(['WARNING: NO UNMASKED VALUES NEAR STATION ',sprintf('%7i',stnID)]);
  xind=1; yind=1;
else
  a=find(rdist==min(rdist));
  xind=xs(a); yind=ys(a);
end
clear dist xs ys x y a ymin ymax xmin xmax

end