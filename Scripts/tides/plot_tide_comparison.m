%======================== plot_tide_comparison.m =========================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  Program plots the Mean Sea Level (MSL) from NOAA CO-OP tide gauge    =%
%=  stations with the predictions from the US East water age model.      =%
%=                                                                       =%
%=========================================================================%

%=========================================================================%
%=                        Load dependent libraries                       =%
%=========================================================================%
addpath(genpath('/he_data/he/zhigang/Matlab_Codes/mexcdf2009/'));
addpath('/home/actodd/MYMATLAB/');

%=========================================================================%
%=                         User-defined indexes                          =%
%=========================================================================%
tStart = datenum(2004,03,01);
tEnd   = datenum(2004,06,01);

%=========================================================================%
%=                   Set Tide file name and load data                    =%
%=========================================================================%
roms.f='/gpfs_share/rhe/actodd/USeast-age/output/useast_his.nc';
tide.f='/home/actodd/ROMS-utils/tides/8729108.txt';

[tide.id]=textread(tide.f,...
    '%7d %*4d-%*2d-%*2d %*2d:%*2d %*7.3f %*7.3f %*d %*d',1,'headerlines',3);
[yy,mm,dd,hh,mi,tide.msl,tide.sigma]=textread(tide.f,...
    '%*7d %4d-%2d-%2d %2d:%2d %7.3f %7.3f %*d %*d','headerlines',3);
tide.date=datenum(yy,mm,dd,hh,mi,mi.*0);

%=========================================================================%
%=                   Set ROMS file name and load data                    =%
%=========================================================================%
roms.f='/gpfs_share/rhe/actodd/USeast-age/output/useast_his.nc';
roms.date=nc_varget(roms.f,'ocean_time',[481],[736]);
roms.lon =nc_varget(roms.f,'lon_rho'               );
roms.lat =nc_varget(roms.f,'lat_rho'               );
%roms.zeta=nc_varget(roms.f,'zeta'      ,[481],[736 1 1]);


%=========================================================================%
%=                      Set list of all station IDs                      =%
%=========================================================================%
List4sID=[8410140;8413320;8418150;8443970;8449130;8447930;8452660;...
          8461490;8465705;8516945;8510560;8531680;8534720;8557380;...
          8570283;8632200;8575512;8636580;8651370;8656483;8658163;...
          8661070;8665530;8670870;8720218;8721604;8722670;8723214;...
          8724580;8723970;8725110;8726384;8726724;8727520;8728690;...
          8729108;8735180;8760922;8764227;8770570;8771341;8772447;...
          8775870;2695540];      
nSta = length(List4sID);

% Tide station information matrix
load('Data/US_TideInfo.mat');

%=========================================================================%
%=                   Download data station by station                    =%
%=========================================================================%
for iSta = 1 : nSta

  %--- Get station index from sID matrix ---%
  inx = find(sID == List4sID(iSta));

  %--- Get the time for last saved record ---%
  saveFile = ['Data/',num2str(ss.ID(inx)),'.dat'];
  if exist(saveFile,'file');
    tmp      = load(saveFile);
    tLatest  = datenum([tmp(end,1:5) 0])+6./60/24;
    clear tmp;
  else
    tLatest  = tStart;
  end
  
  % Download tide observation day by day
  for iT = floor(tLatest) : 1 : floor(tEnd)
    % Time period to donwload in this loop
    tDayGetStart = iT;
    tDayGetEnd   = iT;

    URL = ['http://tidesandcurrents.noaa.gov/api/datagetter?',...
          'product=water_level&',...
          'application=NOS.COOPS.TAC.WL&',...
          'station=',num2str(sID(inx)),'&',...
          'begin_date=',datestr(tDayGetStart,'yyyymmdd'),'&',...
          'end_date=',datestr(tDayGetEnd,'yyyymmdd'),'&',...
          'datum=MSL&',...
          'units=metric&',...
          'time_zone=GMT&',...
          'format=csv'];

    % Grab webpage and save to tmp.dat
    eval(['! wget -O tmp.csv "' URL '"']);
    eval(['! sed -i ''s/' datestr(iT,'yyyy-mm-dd') '/' datestr(iT,'yyyy,mm,dd,') '/g'' tmp.csv']);
    eval(['! sed -i ''s/:/,/g'' tmp.csv']);
    eval(['! sed -i ''s/,p/ /g'' tmp.csv']);
    eval(['! sed -i ''s/,v/ /g'' tmp.csv']);
    
    % Read the data block
    tmp=importdata('tmp.csv',',',1);
    
    % If no valid data then skip
    if ~isfield(tmp,'data') | isempty(tmp.data)
       continue;
    end
    
    % Output array
    Out  = tmp.data(:,1:6);

    % Data quality control
    iBad = abs(Out(:,end)) > 15;
    Out(iBad,end) = nan;

    % Write output
    % Write format define
    saveFormat = '%4d %2d %2d %2d %2d %8.3f\n';
    if exist(saveFile,'file')
       fid = fopen(saveFile,'a');
       fprintf(fid, saveFormat,Out');
       fclose(fid);
    else
       fid = fopen(saveFile,'w');
       fprintf(fid, saveFormat,Out');
       fclose(fid);
    end

    unix(['rm -f tmp2.txt tmp.dat']);
    clear Out tmp yr mth day iBad
  end % End of for iT = floor(tStart) : 1 : floor(tEnd)
  unix(['rm -f tmp.csv']);
end   % for iSta = 1 : nSta       
