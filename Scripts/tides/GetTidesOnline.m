% Get tides from website of NOAA tides&currents: http://tidesandcurrents.noaa.gov/
%
% Written by Zhigang Yao
% North Carolina State University & POL Ocean University of China
% Email: yaozhigang.ouc@gmail.com

% Update history:
%   2013-11-07  : change the download url, to get .csv file, because noaa seems closed 
%                 the previous ascii link.  (Zhigang Yao)

% Station ID list selected
List4sID = [8410140;8413320;8418150;8443970;8449130;8447930;8452660;...
            8461490;8465705;8516945;8510560;8531680;8534720;8557380;...
            8570283;8632200;8575512;8636580;8651370;8656483;8658163;...
            8661070;8665530;8670870;8720218;8721604;8722670;8723214;...
            8724580;8723970;8725110;8726384;8726724;8727520;8728690;...
            8729108;8735180;8760922;8764227;8770570;8771341;8772447;...
            8775870;2695540;];

% Start time to download for the tide observation data
%tStart = datenum([2013 1 1]);
tStart  = floor(now -2);

% End time to download
tEnd   = floor(now-1);

% Save directory for tide observation files
saveDir = './Data/';

% Number of tide station
nSta = length(List4sID);

% Tide station information matrix
load('./Data/US_TideInfo.mat');

for iSta = 1 : nSta

  % Get station index from sID matrix
  inx = find(sID == List4sID(iSta));

  % Get the time for last saved record
  saveFile = fullfile(saveDir,[num2str(ss.ID(inx)) '.dat']);
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

     % Generate corresponding URL
%     URL = ['http://old.tidesandcurrents.noaa.gov/data_listing.shtml?bdate=' datestr(tDayGetStart,'yyyymmdd') '&edate='   datestr(tDayGetEnd,'yyyymmdd')  '&datum=3&unit=0&shift=g&stn=' num2str(sID(inx))  '&type=Tide%20Data&format=View+Data&listing=1'];

      URL = ['http://tidesandcurrents.noaa.gov/api/datagetter?product=water_level&application=NOS.COOPS.TAC.WL&station=' num2str(sID(inx)) '&begin_date=' datestr(tDayGetStart,'yyyymmdd') '&end_date=' datestr(tDayGetEnd,'yyyymmdd') '&datum=MSL&units=metric&time_zone=GMT&format=csv'];

    % Grab webpage and save to tmp.dat
    eval(['! wget -O tmp.csv "' URL '"']);
    eval(['! sed -i ''s/' datestr(iT,'yyyy-mm-dd') '/' datestr(iT,'yyyy,mm,dd,') '/g'' tmp.csv']);
    eval(['! sed -i ''s/:/,/g'' tmp.csv']);
    eval(['! sed -i ''s/Date Time, Water Level, Sigma, O, F, R, L, Quality/2013,11,05,19,36,0.349,0.010,1,0,0,0,p/'' tmp.csv']); 
    eval(['! sed -i ''s/,p/ /g'' tmp.csv']);
    eval(['! sed -i ''s/,v/ /g'' tmp.csv']);
%    tmp = importdata('tmp.csv',' ',1);
   
    % Extract data block, save to tmp2.txt
    % Write a extra line to make sure the column number is right
%    eval(['!echo "8729108 20130224 00 00   0.027   0.147   0.151" > tmp2.txt']);
%    eval(['!grep -ir ^' num2str(ss.ID(inx)) ' tmp.dat | tr : \ >> tmp2.txt']);

    % Read the data block
    tmp = importdata('tmp.csv',',',1);        

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
    saveFormat = '%4d    %2d    %2d    %2d    %2d     %8.3f\n';
    if exist(saveFile,'file')
       fid = fopen(saveFile,'a');
       fprintf(fid, saveFormat, Out');
       fclose(fid);
    else
       fid = fopen(saveFile,'w');
       fprintf(fid, saveFormat, Out');
       fclose(fid);
    end

    unix(['rm -f tmp2.txt tmp.dat']);
    clear Out tmp yr mth day iBad
  end % End of for iT = floor(tStart) : 1 : floor(tEnd)

end   % for iSta = 1 : nSta
