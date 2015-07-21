% This script is to compute the mean age of tracer
% from (a) model history file(s) which have tracer total concentration
% and age content.
%
% Weifeng Gordon Zhang
% July, 2008
%
warning off
hisfile = '/gpfs_share/rhe/SABGOM/FWD/sabgom_avg.nc';
cdlfile = 'sabgom_meanage.cdl';
records = [1:800];   % index of record(s) in hisfile to be converted to ages
ccutoff = 1e-6;      % threshold of minmum tracer concentration 
ageindx = 1;         % index of the release in multiple releases
wrtcoor = 1;         % switch to write coordinate variables into outfile
sec2day = 1;         % convert the age unit from second to day
coorvar = {'spherical','xl','el','theta_s','theta_b','Tcline','hc','s_rho',...
           's_w','Cs_r','Cs_w','h','f','pm','pn','lon_rho','lat_rho','lon_u',...
           'lat_u','lon_v','lat_v','lon_psi','lat_psi','angle','mask_rho',...
           'mask_u','mask_v','mask_psi'};

addpath /home/rhe/
addpath /home/rhe/Matlab_Codes/
startup

outdir  = fileparts(hisfile);
outfile = [outdir '/' cdlfile(1:end-4) '_' num2str(ageindx) '.nc'];
disp(['generating ' outfile ' using ncgen']);
command = ['ncgen -b ' cdlfile];
unix(command);
command = ['mv ' cdlfile(1:end-3) 'nc ' outfile];
unix(command);

if wrtcoor
   command = ['ncks -v '];
   for ic = 1:length(coorvar)
      command = [command coorvar{ic} ','];
   end
   command = [command(1:end-1) ' ' hisfile ' coor_tmp.nc'];
   unix(command);

   command2 = ['ncatted -a  ,global,d,,  coor_tmp.nc'];
   unix(command2);

   command3 = ['ncks -A coor_tmp.nc ' outfile];
   unix(command3);

   command4 = ['rm -f coor_tmp.nc'];
   unix(command4);
end

g     = roms_get_grid(hisfile,hisfile);
mask  = reshape(g.mask_rho,[1 size(g.mask_rho)]);
mask  = repmat(mask,[g.N 1]);

history = ['tracer mean age generated from ' ...
           'history file: ' hisfile ' by Weifeng '...
           ' Zhang using scripts ' which(mfilename) '.m'];

nc_addhist(outfile,history);

for irec = 1:length(records)

   disp(['computing the ages for record :' num2str(records(irec))]);
   % retrieve the data
   rind  = records(irec)-1;
   otime = nc_varget(hisfile,'ocean_time',rind,1);
   nc_varput(outfile,'ocean_time',otime,irec-1,1);

   if ageindx < 6
      varname1 = ['age_0' num2str(2*ageindx-1)];
   else
      varname1 = ['age_' num2str(2*ageindx-1)];
   end
   if ageindx < 5
      varname2 = ['age_0' num2str(2*ageindx)];
   else
      varname2 = ['age_' num2str(2*ageindx)];
   end
   age1  = nc_varget(hisfile,varname1,[rind 0 0 0],[1 -1 -1 -1]);
   age2  = nc_varget(hisfile,varname2,[rind 0 0 0],[1 -1 -1 -1]); 
   meanage = repmat(nan,size(squeeze(age1)));

   % compute the mean age
   ind = find(age1 > ccutoff);
   meanage(ind) = age2(ind) ./ age1(ind);
   if sec2day
      meanage = meanage ./ 86400;
   end

   age1    = reshape(age1,[1 size(age1)]);
   meanage = reshape(meanage,[1 size(meanage)]);

   % write the ages into the netcdf file
   nc_varput(outfile,'age_01',age1,[irec-1 0 0 0],size(age1));
   nc_varput(outfile,'mean_age',meanage,[irec-1 0 0 0],size(age1));
end

timeunit = nc_attget(hisfile,'ocean_time','units');
nc_attput(outfile,'ocean_time','units',timeunit);

disp(['wrote age file: ' outfile]);

