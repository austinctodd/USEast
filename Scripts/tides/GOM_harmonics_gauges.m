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
%=                      Set list of all station IDs                      =%
%=========================================================================%
stnID   =[8724698;8723970;8724580;8725520;8726384;8726724;8727235;...
          8727359;8728130;8728360;8728690;8729108;8729678;8729840;...
          8735180;8742221;8745557;8747437;8760551;8761720;8762075;...
          8764311;8766072;8768094;8771081;8770971;8771510;8772447;...
          8773701;8775237;8775792;8775870;8779750;9500966;8725110;...
          8726520;8727277;8727333;8727306;8727520;8728229;8732828;...
          8741196;8736897;8744117;8747437;8761305;8760943;8761819;...
          8764227;8765251;8771341;8775283;8779770];
stnLonD =[82,81,81,81,82,82,82,82,84,84,84,85,86,87,88,88,89,89,89,89,90,...
          91,92,93,93,94,94,95,96,97,97,97,97,97,81,82,82,82,82,83,84,87,...
          88,88,88,89,89,89,90,91,91,94,97,97];
stnLonM =[55.2,06.3,48.5,52.3,33.9,49.9,38.3,41.5,10.7,30.7,58.9,40.0,...
          51.9,12.7,04.5,40.0,04.9,22.0,08.4,58.1,12.0,23.1,18.3,20.6,...
          38.4,30.8,47.3,18.0,23.3,03.0,14.2,13.0,09.4,47.7,48.4,37.5,...
          41.7,43.4,40.0,01.9,17.4,49.5,32.0,03.5,54.2,19.5,40.4,25.0,...
          02.3,20.4,52.8,43.5,12.2,12.9];
stnLatD =[24,24,24,26,27,27,28,28,30,29,29,30,30,30,30,30,30,30,28,29,...
          29,29,29,29,29,29,29,28,28,27,27,27,26,22,26,27,28,28,28,29,...
          30,30,30,30,30,30,29,28,29,29,29,29,27,26];
stnLatM =[37.9,42.7,33.2,38.8,38.2,58.7,41.5,55.4,04.7,54.9,43.6,09.1,...
          22.6,24.2,15.0,14.3,21.6,16.9,59.4,15.3,06.9,22.3,33.3,45.9,...
          29.9,30.9,17.7,56.0,27.1,49.6,38.0,34.8,04.1,15.7,0.81,45.5,...
          46.3,52.2,51.8,08.1,03.6,25.0,20.4,38.9,24.7,19.5,52.1,55.5,...
          24.1,27.0,42.8,21.5,49.3,03.6];

stnLon=stnLonD+stnLonM./60;
stnLat=stnLatD+stnLatM./60;
      
nSta=length(stnID);

M2amp=[]; M2phase=[]; 
S2amp=[]; S2phase=[];
O1amp=[]; O1phase=[];
K1amp=[]; K1phase=[];
     
for i=1:nSta
  %=======================================================================%
  %=                              Load data                              =%
  %=======================================================================%
  disp(['Reading in data from station ',num2str(stnID(i))]);
  fname=['/home/actodd/ROMS-utils/tides/Data/GOM/',num2str(stnID(i)),'.dat'];
  if exist(fname)
    data=load(fname);

    %=====================================================================%
    %=                               T_tide                              =%
    %=====================================================================%
    consts=['M2  ';'S2  ';'O1  ';'K1  '];
    [tide,xout]=t_tide(data(:,6),'interval',1/10,...
                                 'start time',[2004,1,4,0],...
                                 'latitude',stnLat(i));
    for k=1:size(tide.name,1)
      if (strcmp(tide.name(k,:),consts(1,:))==1)
        M2amp(  i)=tide.tidecon(k,1);
        M2phase(i)=tide.tidecon(k,3);
      elseif (strcmp(tide.name(k,:),consts(2,:))==1)
        S2amp(  i)=tide.tidecon(k,1);
        S2phase(i)=tide.tidecon(k,3);
      elseif (strcmp(tide.name(k,:),consts(3,:))==1)
        O1amp(  i)=tide.tidecon(k,1);
        O1phase(i)=tide.tidecon(k,3);
      elseif (strcmp(tide.name(k,:),consts(4,:))==1)
        K1amp(  i)=tide.tidecon(k,1);
        K1phase(i)=tide.tidecon(k,3);
      end
    end
  else
    disp('FILE DOES NOT EXIST');
  end
end
  
save('/home/actodd/ROMS-utils/tides/Data/GOM/GOM_harmonics.mat',...
     'M2phase','M2amp','S2phase','S2amp',...
     'O1phase','O1amp','K1phase','K1amp','-MAT');
 
 


