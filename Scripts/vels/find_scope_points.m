%=============================================================================%
%=                                                                           =%
%=  Created by Austin C. Todd, NCSU, 2015                                    =%
%=                                                                           =%
%=  Program is designed to find the points on the US East Grid where at the  =%
%=  outer boundary of the shelf residence time scope mask.                   =%
%=                                                                           =%
%=============================================================================%

%--- Startup and setting of the path ---%
addpath('~/MYMATLAB/');
addpath(genpath('~/MYMATLAB/ROMS-matlab/'));
format long g; format compact;

%=============================================================================%
%=                                                                           =%
%=       !!!!!!!!!!!        USER-DEFINED VARIABLES      !!!!!!!!!!           =%
%=                                                                           =%
%=============================================================================%
ATL_file='/Volumes/Black_box/Data/USeast/data/grd/grid_ATLscope.nc';
GOM_file='/Volumes/Black_box/Data/USeast/data/grd/grid_GOM_shelf_scope.nc';
out_file='/Users/todd/Documents/Work/Projects/USeast-age/shelf_points.txt';

%=============================================================================%
%=                  Setup some grid and preliminary arrays                   =%
%=============================================================================%
ncid=netcdf.open(ATL_file,'nowrite');
    h     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'         ));
    lon   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'   ));
    lat   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'   ));
    hmask = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'));
    Ascope= netcdf.getVar(ncid,netcdf.inqVarID(ncid,'scope_rho'  ));
netcdf.close(ncid)

ncid=netcdf.open(GOM_file,'nowrite');
    Gscope= netcdf.getVar(ncid,netcdf.inqVarID(ncid,'scope_rho'  ));
netcdf.close(ncid)

Sscope=Ascope+Gscope;

%=============================================================================%
%=                Find scope points at outer shelf boundary                  =%
%=============================================================================%
[c,hc]=contour(Sscope,[1 1],'k');
hh    = []; hs    = []; hi    = []; hj    = [];
hh100 = []; hs100 = []; hi100 = []; hj100 = [];
count = 1; count100 = 1;
for i=1:length(c)
    if (c(1,i)>5 && c(2,i)>5 && c(1,i)<482 && c(2,i)<402)
        if (h(c(2,i),c(1,i))>100)
            hh100(count100)=h(     c(2,i),c(1,i));
            hs100(count100)=Sscope(c(2,i),c(1,i));
            hi100(count100)=c(1,i); 
            hj100(count100)=c(2,i);
            count100=count100+1;
        end
        hh(count)=h(     c(2,i),c(1,i));
        hs(count)=Sscope(c(2,i),c(1,i));
        hi(count)=c(1,i); 
        hj(count)=c(2,i);
        count=count+1;
    end
end

% Cleanup those points that are still on land boundary (or Bermuda)
hh100=hh100(1:722);
hs100=hs100(1:722);
hi100=hi100(1:722);
hj100=hj100(1:722);

% Write indexes to file
fid=fopen(out_file,'w')
for i=1:722
    fprintf(fid,'%03i %03i\n',[hj100(i)-1,hi100(i)-1]);
end
fclose(fid);

return;
