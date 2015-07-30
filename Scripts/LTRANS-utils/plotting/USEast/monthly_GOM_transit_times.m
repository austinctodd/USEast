%=========================================================================%
%
% PROGRAM: GOM_transit_times.m
%
% PURPOSE: Program reads in LTRANS output and calculates the length of time
%          that each particle spends within the specified residence time
%          domain. Transit/Exposure times are calculated for each month.
%
% AUTHOR: Austin C. Todd (NCSU), 9 March 2015
%
%=========================================================================%

clear all
%=========================================================================%
% Set file directories
%=========================================================================%
LTRANS_dir = '/Volumes/Black_box/Data/LTRANS/output/Mississippi/';
grid_file  =['/Volumes/Black_box/Data/USeast/Data/grd/',....
             'grid_GOM_shelf_scope.nc'];

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

day_count=[31,28,31,30,31,30,31,31,30,31,30,31];

for exps=1:5

  clear plon plat pdepth pcol pdob plonind platind pscope   

  %========================================================================
  % Create file names for particle file and grid file
  %========================================================================
  particle_file =[LTRANS_dir,'AR.nc'];%',sprintf('%1i',exps),'.nc'];
  
  %========================================================================
  %  Open particle file and read in data
  %========================================================================
  disp(['Reading data from LTRANS file']);
  ncid=netcdf.open(particle_file,'nowrite');
    plon   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
    plat   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
    pdepth = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'depth'));
    pcol   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'color'));
    pdob   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dob'));
  netcdf.close(ncid);
  pdob=floor(pdob./86400);
  plonind=floor((plon-lon(1,1))*10.2790697675)+1;

  %========================================================================
  % Find longitude indexes from particle longitudes
  %========================================================================
  disp(['Calculating global particle stats...']);
  for i=1:size(plon,1)
    for j=pdob(i)+1:size(plon,2)
      platind(i,j)=max(find(lat(1,:)<=plat(i,j)));
      if (pcol(i,j)>-1)
        pscope( i,j)=scope(plonind(i,j),platind(i,j));
      else
        pscope( i,j)=0;
      end
    end
  end

  %======================================================================
  % Find longitude indexes from particle longitudes
  %======================================================================
  pdens      =zeros(12,402,482); 
  trans_times=zeros(12,402,482);
  trans_parts=zeros(12,402,482);
  expos_times=zeros(12,402,482);
return;
  for i=1:size(plon,1)      
    for month=1:12
      disp(['Calculating particle stats for Month ',smn]);
      smn=sprintf('%02i',month);
      inds=sum([0 day_count(month)]);
      
      pind=[inds inds+365 inds+730];
      
      
      for j=pdob(i)+1:size(plon,2)
        if (pcol(i,j)>-1)
          pdens(plonind(i,j),platind(i,j))=pdens(plonind(i,j),platind(i,j))+1;
        end
      end
      ind=min(find(pscope(i,pdob(i)+1:end)<1));
      if (length(ind)<1)   
        if pscope(i,pdob(i)+2)==1
          ind2=size(plon,2);
        else
          ind2=pdob(i);
        end
      else
        ind2=pdob(i)+ind;
      end
      part_ttime(i)=ind2-pdob(i)+1;
      part_expos(i)=length(find(pscope(i,pdob(i)+1:end)==1));
      for j=pdob(i)+1:ind2        
        %==================================================================
        % Transit time before leaving shelf
        %==================================================================
        trans_times(month,plonind(i,j),platind(i,j))=...
        trans_times(month,plonind(i,j),platind(i,j))+(ind2-pdob(i)+1);
        
        %==================================================================
        % Exposure time (total length particle stays on shelf)
        %==================================================================
        expos_times(month,plonind(i,j),platind(i,j))=...
        expos_times(month,plonind(i,j),platind(i,j))+...
                    length(find(pscope(i,pdob(i)+1:end)==1));
                
        %==================================================================
        % Add number of partilces through this location
        %==================================================================
        trans_parts(month,plonind(i,j),platind(i,j))=...
        trans_parts(month,plonind(i,j),platind(i,j)) +1;
      end
    end

  for month=1:12
    %======================================================================
    % Output gridded mean transport times to file
    %======================================================================
    output_file   =[LTRANS_dir,'MS',sprintf('%1i',exps),'_',smn,'.txt'];
    fid=fopen(output_file,'w');
    for i=1:200
      for j=100:260
        fprintf(fid,'%03i %03i %12.5f %12.5f %05i\n',...
                [i j trans_times(i,j) expos_times(i,j) trans_parts(i,j)]);
      end
    end
    fclose(fid);
    

  end
end
