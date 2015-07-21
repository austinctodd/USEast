% This script is to locate the spill rig in SABGOM grid and compute
% the oil concentration increase rate at the release sites (point
% release or column release)
%
% By Weifeng Zhang
% Aug. 2010

riglat  = 28.736628;
riglon  = -88.365997;
flowrate= 9857;    % official estimate of the flow rate on Aug. 2 (m^3 / day)
prelease= 0;       % switch for point release or column release
                   % 1 -- point release;  0 -- column release
depth   = [-1200 -500 -5];   % only valid for point release

hisfile = '/gpfs_share/rhe/SABGOM/FWD/sabgom_his.nc';

g       = roms_get_grid(hisfile,hisfile);

[j,i]   = closest(g.lon_rho,g.lat_rho,[riglon riglat]);

if prelease
   % for point release
   kk      = repmat(nan,size(depth));
   conrate = repmat(nan,size(depth));

   for id  = 1:length(depth)
      tind = find(g.z_w(:,j,i) < depth(id));
      kk(id) = tind(end);

      cellhght= g.z_w(kk(id)+1,j,i)-g.z_w(kk(id),j,i);
      cellvol = 1./g.pm(j,i)./g.pn(j,i).*cellhght;
      trpt    = flowrate./86400;
      conrate(id) = trpt./cellvol;
   end
   disp(' ');
   disp(['Layer indices of the release points are:']);
   disp(['  ' num2str(kk)]);
   disp(' ');
   disp(['Oil concentration increase rates at the release points are:']);
   disp(['  ' num2str(conrate)]);

else
   wdpth  = g.h(j,i);
   colvol = 1./g.pm(j,i)./g.pn(j,i).*wdpth;
   trpt   = flowrate./86400;
   conrate= trpt./colvol; 

   disp(' ');
   disp(['Oil concentration increase rate for the whole water column is:']);
   disp(['  ' num2str(conrate)]);

end

