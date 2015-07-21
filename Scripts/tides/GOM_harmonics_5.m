%======================== USeast-tidal-harmonics =========================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program is designed to read in results from a tidal-forced US   =%
%=  East Coast ROMS simulation and perform tidal harmonic analysis on    =%
%=  the simulated sea surface elevation.  The program is inherently      =%
%=  dependent on the t_tide routines located in the directory:           =%
%=               /home/actodd/MYMATLAB/ROMS-matlab/t_tide/               =%
%=                                                                       =%
%=========================================================================%

%--- Add various libraries and paths ---%
addpath(genpath('/he_data/he/zhigang/Matlab_Codes/mexcdf2009'));
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));

%=========================================================================%
%=                               Load data                               =%
%=========================================================================%
disp('Reading in ROMS output');
fname='/gpfs_share/rhe/actodd/USeast-tide/output/LTP_OBC/newcode_orig_TDATE/useast_his.nc';
zeta=nc_varget(fname,'zeta');
lat =nc_varget(fname,'lat_rho');
mask=nc_varget(fname,'mask_rho');

%=========================================================================%
%=                                 T_tide                                =%
%=========================================================================%
consts=['M2  ';'S2  ';'O1  ';'K1  '];
for i=81:263
  for j=10:186
    disp(['i = ',sprintf('%03i',i),' , j = ',sprintf('%03i',j)]);
    if (mask(i,j)==1)
      [tide,xout]=t_tide(zeta(97:745,i,j),...
                         'interval',1,...
                         'start time',[2004,1,5,0],...
                         'latitude',lat(i,j));
      for k=1:size(tide.name,1)
        if (strcmp(tide.name(k,:),consts(1,:))==1)
          M2amp(  i-80,j-9)=tide.tidecon(k,1);
          M2phase(i-80,j-9)=tide.tidecon(k,3);
        elseif (strcmp(tide.name(k,:),consts(2,:))==1)
          S2amp(  i-80,j-9)=tide.tidecon(k,1);
          S2phase(i-80,j-9)=tide.tidecon(k,3);
        elseif (strcmp(tide.name(k,:),consts(3,:))==1)
          O1amp(  i-80,j-9)=tide.tidecon(k,1);
          O1phase(i-80,j-9)=tide.tidecon(k,3);
        elseif (strcmp(tide.name(k,:),consts(4,:))==1)
          K1amp(  i-80,j-9)=tide.tidecon(k,1);
          K1phase(i-80,j-9)=tide.tidecon(k,3);
        end
     end
    else
      M2phase(i-80,j-9)=nan; M2amp(i-80,j-9)=nan;
      S2phase(i-80,j-9)=nan; S2amp(i-80,j-9)=nan;
      O1phase(i-80,j-9)=nan; O1amp(i-80,j-9)=nan;
      K1phase(i-80,j-9)=nan; K1amp(i-80,j-9)=nan;
    end
  end
end
save(['/gpfs_share/rhe/actodd/USeast-tide/output/LTP_OBC/newcode_orig_TDATE/',...
      'GOM_harmonic_analysis.mat'],'M2phase','M2amp','S2phase','S2amp',...
      'O1phase','O1amp','K1phase','K1amp','-MAT');
 
 


