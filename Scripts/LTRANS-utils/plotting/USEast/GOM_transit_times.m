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
for exps=3:5
clear plon plat pdepth pcol pdob plonind platind pscope   
    
%=========================================================================%
% Set file directories
%=========================================================================%
LTRANS_dir ='/Volumes/Black_box/Data/LTRANS/output/Mississippi/';
ROMS_dir   ='/Volumes/Black_box/Data/USeast/Data/grd/';

%=========================================================================%
% Create file names for particle file and grid file
%=========================================================================%
particle_file =[LTRANS_dir,'AR.nc'];%',sprintf('%1i',exps),'.nc'];
grid_file     =[ROMS_dir,'grid_GOM_shelf_scope.nc'];
output_file   =[LTRANS_dir,'AR.txt'];%',sprintf('%1i',exps),'.txt'];

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

disp(['Looping through 29230 particles...']);
for i=1:size(plon,1)
  for j=pdob(i)+1:size(plon,2)
    platind(i,j)=max(find(lat(1,:)<=plat(i,j)));
    if (pcol(i,j)>-1)
%        plonind(i,j)=round(interp1(lon(:,1),[1:402],plon(i,j),'nearest'));
%        platind(i,j)=round(interp1(lat(1,:),[1:482],plat(i,j),'nearest'));
%       pscope( i,j)=ceil(interp2(lon,lat,scope,plon(i,j),plat(i,j),'nearest'));
        pscope( i,j)=scope(plonind(i,j),platind(i,j));
        pdens(plonind(i,j),platind(i,j))=pdens(plonind(i,j),platind(i,j))+1;
    else
        pscope( i,j)=0;
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
      trans_times(plonind(i,j),platind(i,j))=...
          trans_times(plonind(i,j),platind(i,j))+(ind2-pdob(i)+1);
      expos_times(plonind(i,j),platind(i,j))=...
          expos_times(plonind(i,j),platind(i,j))+...
          length(find(pscope(i,pdob(i)+1:end)==1));
      trans_parts(plonind(i,j),platind(i,j))=...
          trans_parts(plonind(i,j),platind(i,j))+1;
  end
end

%=========================================================================%
% Output gridded mean transport times to file
%=========================================================================%
fid=fopen(output_file,'w');
for i=1:200
  for j=100:260
      fprintf(fid,'%03i %03i %12.5f %12.5f %05i\n',...
              [i j trans_times(i,j) expos_times(i,j) trans_parts(i,j)]);
  end
end
fclose(fid);
end

return;
%=========================================================================%
% Make crude animation of just particles remaining on shelf
%=========================================================================%
figure(2)
plotdir='/Volumes/Black_box/Data/PLOTS/LTRANS/Mississippi/frames/';
for i=1:365
  aa = find(pcol(:,i)>-1);
  alon   = plon(aa,i);
  alat   = plat(aa,i);
  ascope = pscope(aa,i);
  aage   = i-pdob(aa);
  atrans = part_ttime(aa);
  aexpos = part_etime(aa);
  
  aa = find(ascope>0);
  alon   = alon(aa);
  alat   = alat(aa);
  aage   = aage(aa);
  atrans = atrans(aa);
  aexpos = aexpos(aa);
  
  clf
  contourf(lon,lat,1-mask,[0 1],'k')
  caxis([0 3]); colormap(1-gray)
  axis xy; axis image
  axis([-98 -80 18 31])
  hold on;
  contour(lon,lat,scope,[0 0],'k');
  freezeColors
  colormap(jet)
  scatter(alon,alat,5,atrans,'filled')
  title(['Day ',sprintf('%04i',i)]);
  caxis([0 365.25])
  colorbar;

  print('-dpng','-r150','-painters',[plotdir,'frame_',sprintf('%03i',i),'.png']); 
  pause(0.25);
end
