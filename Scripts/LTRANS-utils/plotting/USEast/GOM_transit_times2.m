%=========================================================================%
%
% PROGRAM: GOM_transit_times.m
%
% PURPOSE: Program reads in LTRANS output and calculates the length of time
%          that each particle spends within the specified residence time
%          domain. 
%
% AUTHOR: Austin C. Todd (NCSU), 9 March 2015
%
%=========================================================================%

clear all
for exps=1:1
clear plon plat pdepth pcol pdob plonind platind pscope   
clear part_ttime part_etime

%=========================================================================%
% Set file directories
%=========================================================================%
LTRANS_dir ='/Volumes/Black_box/Data/LTRANS/output/Mississippi/';
ROMS_dir   ='/Volumes/Black_box/Data/USeast/Data/grd/';

%=========================================================================%
% Create file names for particle file and grid file
%=========================================================================%
particle_file =[LTRANS_dir,'MS',sprintf('%1i',exps),'.nc'];
grid_file     =[ROMS_dir,'grid_GOM_shelf_scope.nc'];
output_file   =[LTRANS_dir,'MS',sprintf('%1i',exps),'_v2.txt'];

%=========================================================================%
%  Open particle file and read in data
%=========================================================================%
disp(['Reading data from LTRANS file']);
ncid=netcdf.open(particle_file,'nowrite');
  plon   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
  plat   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
  pdepth = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'depth'));
  pcol   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'color'));
  pdob   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dob'));
netcdf.close(ncid);
pdob=floor(pdob./86400);

%=========================================================================%
% Open grid file and read lat/lon/scope
%=========================================================================%
disp(['Reading data from ROMS file']);
ncid=netcdf.open(grid_file,'nowrite');
  lon   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'));
  lat   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'));
  scope = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'scope_rho'));
  mask  = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'));
netcdf.close(ncid);

%=========================================================================%
% Find longitude indexes from particle longitudes
%=========================================================================%
plonind=floor((plon-lon(1,1))*10.2790697675)+1;
pdens=zeros(402,482); 
trans_times=zeros(402,482);
trans_parts=zeros(402,482);
expos_times=zeros(402,482);
expos_parts=zeros(402,482);

disp(['Looping through 29230 particles...']);
for i=1:size(plon,1)
  for j=pdob(i)+1:size(plon,2)
    platind(i,j)=max(find(lat(1,:)<=plat(i,j)));
    if (pcol(i,j)>-1)
        pscope( i,j)=scope(plonind(i,j),platind(i,j));
    else
        pscope( i,j)=0;
    end
  end
  
  %========================================================================
  % Find first time particle exits the shelf
  %========================================================================
  ind1=min(find(pscope(i,pdob(i)+1:end)<1));
  ind2=max(find(pscope(i,pdob(i)+1:end)>0));
   
  if (length(ind1)<1)   %<-- Particle never leaves the shelf
    if pscope(i,pdob(i)+2)==1
      t_ind=size(plon,2);
      e_ind=size(plon,2);
      part_ttime(i)=size(plon,2)-pdob(i)+1;
      part_etime(i)=size(plon,2)-pdob(i)+1;
    else
      t_ind=pdob(i);
      e_ind=pdob(i);
      part_ttime(i)=1;
      part_etime(i)=1;
    end
  else
    part_ttime(i)=ind1;
    t_ind=ind1;
    if (ind2<ind1) %<-- Particle never re-enters the shelf
      e_ind=pdob(i)+ind1;
      part_etime(i)=ind1;
    else           %<-- Particle re-enters the shelf at some point
      e_ind=pdob(i)+ind2;
      part_etime(i)=length(find(pscope(i,pdob(i)+1:end)==1));
    end
  end
  
  %========================================================================
  % Map the transit & residence times to the model grid
  %========================================================================
  for j=pdob(i)+1:e_ind
    if (pscope(i,j)==1)
      if j<=t_ind
        trans_times(plonind(i,j),platind(i,j))=...
        trans_times(plonind(i,j),platind(i,j))+part_ttime(i);
        trans_parts(plonind(i,j),platind(i,j))=...
        trans_parts(plonind(i,j),platind(i,j))+1;
        expos_times(plonind(i,j),platind(i,j))=...
        expos_times(plonind(i,j),platind(i,j))+part_etime(i);
        expos_parts(plonind(i,j),platind(i,j))=...
        expos_parts(plonind(i,j),platind(i,j))+1;
      else
        expos_times(plonind(i,j),platind(i,j))=...
        expos_times(plonind(i,j),platind(i,j))+part_etime(i);
        expos_parts(plonind(i,j),platind(i,j))=...
        expos_parts(plonind(i,j),platind(i,j))+1;
      end
    end
  end
end
return;
%=========================================================================%
% Output gridded mean transport times to file
%=========================================================================%
%output_file   =[LTRANS_dir,'AR_v2.txt'];%',sprintf('%1i',exps),'.txt'];
fid=fopen(output_file,'w');
for i=1:200
  for j=100:260
      fprintf(fid,'%03i %03i %12.5f %05i %12.5f %05i\n',...
              [i,j,trans_times(i,j),trans_parts(i,j),...
                   expos_times(i,j),expos_parts(i,j)]);
  end
end
fclose(fid);

end

