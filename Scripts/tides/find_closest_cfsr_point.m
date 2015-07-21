function [xind,yind]=find_closest_cfsr_point(cfsr,ss,stnID);
%FIND_CLOSEST_CFSR_POINT
%=========================================================================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  Program finds the ROMS grid lat/lon pairs that are closest to a tide =%
%=  gauge station given by ID stnID                                      =%
%=                                                                       =%
%=  INPUT VARIABLES:                                                     =%
%=    cfsr  ::  structured array with variables cfsr.lon and cfsr.lat    =%
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
xmax=min(find(cfsr.lon>=ss.lon(f_ind)));
xmin=max(find(cfsr.lon<=ss.lon(f_ind)));
ymin=max(find(cfsr.lat<=ss.lat(f_ind)));
ymax=min(find(cfsr.lat>=ss.lat(f_ind)));

count=1;
for x=xmin:xmax
  for y=ymin:ymax
    dist(count)=sw_dist([ss.lat(f_ind) cfsr.lat(y)],...
                        [ss.lon(f_ind) cfsr.lon(x)]);
    xs(count)=x; ys(count)=y;
    count=count+1;  
  end
end
a=find(dist==min(dist));
xind=xs(a); yind=ys(a);
clear dist xs ys x y a ymin ymax xmin xmax

end