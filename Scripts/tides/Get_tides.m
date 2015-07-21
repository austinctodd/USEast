%============================= Get_tides.m ===============================%
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
tStart = datenum(2004,01,01);
tEnd   = datenum(2006,01,01);

%=========================================================================%
%=                      Set list of all station IDs                      =%
%=========================================================================%
List4sID=[8410140;8413320;8418150;8443970;8449130;8447930;8452660;...
          8461490;8465705;8516945;8510560;8531680;8534720;8557380;...
          8570283;8632200;8575512;8636580;8651370;8656483;8661070;...
          8665530;8670870;8720218;8721604;8723214;8760922;8770570;...
	  2695540;8724698;8723970;8724580;8725520;8726384;8726724;...
	  8727235;8727359;8728130;8728360;8728690;8729108;8729678;...
	  8729840;8735180;8742221;8745557;8747437;8760551;8761720;...
	  8762075;8764311;8766072;8768094;8771081;8770971;8771510;...
	  8772447;8773701;8775237;8775792;8775870;8779750;9500966;...
	  8725110;8726520;8727277;8727333;8727306;8727520;8728229;...
	  8732828;8741196;8736897;8744117;8747437;8761305;8760943;...
	  8761819;8764227;8765251;8771341;8775283;8779770];
nSta = length(List4sID);

% Tide station information matrix
load('Data/US_TideInfo.mat');

%=========================================================================%
%=                   Download data station by station                    =%
%=========================================================================%
for iSta = 1 : nSta

  disp(['Station ',num2str(List4sID(iSta))])
  
  %--- Get the time for last saved record ---%
  saveFile = ['/he_data/he/actodd/DATA/tides/',num2str(List4sID(iSta)),'.dat'];
  if exist(saveFile,'file');
    tmp      = load(saveFile);
    tLatest  = datenum([tmp(end,1:5) 0])+6./60/24;
    clear tmp;
  else
    tLatest  = tStart;
  end
  
  % Download tide observation day by day
  for iT = floor(tLatest):floor(tEnd)
    % Time period to donwload in this loop
    tDayGetStart = iT;
    tDayGetEnd   = iT;

    URL = ['http://tidesandcurrents.noaa.gov/api/datagetter?',...
          'product=water_level&',...
          'application=NOS.COOPS.TAC.WL&',...
          'station=',num2str(List4sID(iSta)),'&',...
          'begin_date=',datestr(tDayGetStart,'yyyymmdd'),'&',...
          'end_date=',datestr(tDayGetEnd,'yyyymmdd'),'&',...
          'datum=MSL&',...
          'units=metric&',...
          'time_zone=GMT&',...
          'format=csv'];

    % Grab webpage and save to tmp.dat
    eval(['! wget -O tmp.csv "' URL '"']);
    eval(['! sed -i ''s/' datestr(iT,'yyyy-mm-dd'),'/',...
                          datestr(iT,'yyyy,mm,dd,') '/g'' tmp.csv']);
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

 
end   % for iSta = 1 : nSta       
