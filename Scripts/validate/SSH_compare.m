%============================ SSH_compare.m ==============================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program is designed to read in results from a the US East Coast =%
%=  water age ROMS simulation and calculate the transport through the FL =%
%=  Straits at 27N, and compare to observed transport values. This       =%
%=  program requires use of the ROMS grid stretching function set_depth  =%
%=                                                                       =%
%=========================================================================%

%--- Add various libraries and paths ---%
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));

%=========================================================================%
%=                       Load AVISO data for 2004                         =%
%=========================================================================%
disp('Reading in ROMS output');
load('/he_data/he/actodd/DATA/USeast/ssh/aviso_ssh_2004-2006.mat');
lon=0:0.25:



%=========================================================================%
%=                       Load ROMS data for 2005                         =%
%=========================================================================%
pdir='/gpfs_share/rhe/actodd/USeast-age/output/OBC/2005/';

for i=367:629
 disp(['2005 :: day ',sprintf('%3i',i)]);
 if i==367
   fname=[pdir,'useast_his_',sprintf('%04i',i),'.nc'];
   ncid=netcdf.open(fname,'nowrite');
   zeta(:,i)=double(squeeze(netcdf.getVar(ncid,netcdf.inqVarID(ncid,...
                                      'zeta'),[195 208   0],[11 1    1])));
   v( :,:,i)=double(squeeze(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                                    [195 207 0 0],[11 1 36 1])));  
   netcdf.close(ncid);
 else
   fname=[pdir,'useast_his_',sprintf('%04i',i-1),'.nc'];
   ncid=netcdf.open(fname,'nowrite');
   zeta(:,i)=double(squeeze(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'),...
                                     [195 208   7],[11 1    1])));
   v( :,:,i)=double(squeeze(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                                     [195 207 0 7],[11 1 36 1])));
   netcdf.close(ncid);
 end
 zw=squeeze(set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,...
                  36,5,h,zeta(:,i),0));
 zw=zw(:,2:end)-zw(:,1:end-1); 
 transport(:,i)=nansum(zw.*squeeze(v(:,:,i)),2).*dist;
end

pdir='/gpfs_share/rhe/actodd/USeast-age/output/';
for i=640:988
 disp(['2005 :: day ',sprintf('%3i',i)]);
 fname=[pdir,'useast_his_',sprintf('%04i',i-1),'.nc'];
 ncid=netcdf.open(fname,'nowrite');
   zeta(:,i)=double(squeeze(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'),...
                                     [195 208   7],[11 1    1])));
   v( :,:,i)=double(squeeze(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                                     [195 207 0 7],[11 1 36 1])));
 netcdf.close(ncid);
 zw=squeeze(set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,...
                  36,5,h,zeta(:,i),0));
 zw=zw(:,2:end)-zw(:,1:end-1); 
 transport(:,i)=nansum(zw.*squeeze(v(:,:,i)),2).*dist;
end
%=========================================================================%
%=                            Load ROMS data                             =%
%=========================================================================%
disp('Reading in FL Current data');
[FC.yr,FC.mn,FC.dy,FC.transport,FC.flag]=textread(['data/',...
    'FC_cable_transport_2004.dat'],'%4d %2d %2d %f %1d','headerlines',21);
[FC.yr(end+1:end+365),FC.mn(end+1:end+365),FC.dy(end+1:end+365),...
 FC.transport(end+1:end+365),FC.flag(end+1:end+365)]=textread(['data/',...
    'FC_cable_transport_2005.dat'],'%4d %2d %2d %f %1d','headerlines',21);
[FC.yr(end+1:end+365),FC.mn(end+1:end+365),FC.dy(end+1:end+365),...
 FC.transport(end+1:end+365),FC.flag(end+1:end+365)]=textread(['data/',...
    'FC_cable_transport_2006.dat'],'%4d %2d %2d %f %1d','headerlines',33);


return;

%=========================================================================%
%=                                 T_tide                                =%
%=========================================================================%
