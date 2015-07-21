% Plot taylor diagrams for each station


%add project information
run ./loadUSeastInfo.m;

%
staMatFile = '/raid0/zyao/Project/USeast/Scripts/Proj_JOO/station.mat';

% Load mat file
load(staMatFile);
tStart = tROMS(1);
tEnd   = tROMS(end);
tTarget = tStart : 3./24 : tEnd;
xTick   = tStart :   5   : tEnd;


%
%ha = tight_subplot(1,3,[.05 .05],[.1 .05],[.05 .05]);
set(gcf, 'Color', 'w');

obsID = {  ...
    '8410140'; ...
    '8413320'; ...
    '8418150'; ...
    '8443970'; ...
    '8449130'; ...
    '8447930'; ...
    '8452660'; ...
    '8461490'; ...
    '8465705'; ...
    '8516945'; ...
    '8510560'; ...
    '8531680'; ...
    '8534720'; ...
    '8557380'; ...
    '8570283'; ...
    '8632200'; ...
    '8575512'; ...
    '8636580'; ...
    '8651370'; ...
    '8656483'; ...
    '8658163'; ...
    '8661070'; ...
    '8665530'; ...
    '8670870'; ...
    '8720218'; ...
    '8721604'; ...
    '8722670'; ...
    '8723214'; ...
    '8724580'; ...
    '8723970'; ...
    '8725110'; ...
    '8726384'; ...
    '8726724'; ...
    '8727520'; ...
    '8728690'; ...
    '8729108'; ...
    '8735180'; ...
    '8760922'; ...
    '8764227'; ...
    '8770570'; ...
    '8771341'; ...
    '8772447'; ...
    '8775870'; ...
    '2695540'; ...
    };

% =============================================================================%

% Number of observations
nObs = length(obsID);

% NOS Tide observation directory
obsDir  = '/raid0/zyao/Datasets/Tides/Data/';

% Tidal observation informat matrix
load('/raid0/zyao/Datasets/Tides/Scripts/Info.mat');


for iObs = 1 : nObs
    clear obsTide obsTime
    
    % Read observation
    obsFile = fullfile(obsDir,[obsID{iObs} '.dat']);
    obsTide = load(obsFile);
    obsTime = datenum([obsTide(:,1:5) zeros(size(obsTide,1),1)]);
    
    
    % Model results
    modTime = tROMS;
    modTide = staROMS(:,iObs);
    
    
    % Find good observation point lying the plot window
    iKeep = find(obsTime >= tStart & obsTime <= tEnd & ~isnan(obsTide(:,end)));
    obsNew = obsTide(iKeep,end);
    % Map the modeled result onto the same point
    modNew = interp1(modTime, modTide, obsTime(iKeep));
    
    %  Error=Error1(T,P,Comment);
    tmp =  Error1(obsNew-nanmean(obsNew), modNew-nanmean(modNew),'');
    stdTide(iObs,1) = tmp.StdP./tmp.StdT;
    rmsTide(iObs,1) = tmp.RMSE;
    corTide(iObs,1) = tmp.Cor;
    
    
end
% Remove two bad points
iGo = corTide < 0.5;
stdTide(iGo) = [];
rmsTide(iGo) = [];
corTide(iGo) = [];
nStn4Tide    = length(rmsTide);

if 0
    stdTide = [1; stdTide(:)];
    rmsTide = [0; rmsTide(:)];
    corTide = [1; corTide(:)];
end

%axes(ha(1));
%subplot(1,3,1);
%[hp,ht,axl] = taylordiag_USeast(stdTide,rmsTide,corTide,'tickRMS',[0.5:0.5:1.5]);

% =============================================================================%
%                   Start Weather Plot for each station                          %
% =============================================================================%

obsID = {  ...
    41002;...
    41008;...
    41012;...
    41013;...
    41046;...
    41047;...
    41048;...
    42001;...
    42003;...
    42012;...
    42020;...
    42036;...
    42056;...
    42058;...
    44014;...
    44025;...
    44037;...
    44056;...
    44065;...
    99999;...
    };

% NOS Tide observation directory
obsDir  = '/raid0/zyao/Datasets/Meteorological/Data/';
nObs = length(obsID);

clear stdPair corPair rmsPair
for  iObs = 1: nObs;
    
    sIDName = num2str(obsID{iObs});
    if obsID{iObs} == 9999
        sIDName = 'LOPL1';
    end
    
    % Observation file
    obsFile = fullfile(obsDir,[num2str(obsID{iObs}) '.dat']);
    
    obsATM = load(obsFile);
    obsTime = datenum([obsATM(:,1:5) zeros(size(obsATM,1),1)]); otInx= 1:length(obsTime);
    tObs    = obsTime(otInx);
    obsUwind= -1*obsATM(otInx,7).*sind(obsATM(otInx,6));
    obsVwind= -1*obsATM(otInx,7).*cosd(obsATM(otInx,6));
    obsHwave= obsATM(otInx,9);
    obsDwave= obsATM(otInx,10);
    obsPair = obsATM(otInx,13);
    
    % Modeled
    tmpData  = squeeze(staWRF(iObs,:,:));
    modTime  = tWRF;
    modTair  = tmpData(:,1);
    modPair  = tmpData(:,2);
    modUwind = tmpData(:,3);
    modVwind = tmpData(:,4);
    modHwave = tmpData(:,5);
    
    % Find good observation point lying the plot window
    iKeep = find(obsTime >= tStart & obsTime <= tEnd & ~isnan(obsPair(:,end)));
    obsNew = obsPair(iKeep,end);
    % Map the modeled result onto the same point
    modNew = interp1(modTime, modPair, obsTime(iKeep));
    
    %  Error=Error1(T,P,Comment);
    tmp =  Error1(obsNew, modNew,'');
    stdPair(iObs,1) = tmp.StdP./tmp.StdT;
    rmsPair(iObs,1) = tmp.RMSE;
    corPair(iObs,1) = tmp.Cor;
    
    % Find good observation point lying the plot window
    iKeep = find(obsTime >= tStart & obsTime <= tEnd & ~isnan(obsHwave(:,end)));
    obsNew = obsHwave(iKeep,end);
    % Map the modeled result onto the same point
    modNew = interp1(modTime, modHwave, obsTime(iKeep));
    
    %  Error=Error1(T,P,Comment);
    tmp =  Error1(obsNew, modNew,'');
    stdHwave(iObs,1) = tmp.StdP./tmp.StdT;
    rmsHwave(iObs,1) = tmp.RMSE;
    corHwave(iObs,1) = tmp.Cor;
    
    
    % Find good observation point lying the plot window
    iKeep = find(obsTime >= tStart & obsTime <= tEnd & ~isnan(obsUwind(:,end)));
    obsNew = obsUwind(iKeep,end);
    % Map the modeled result onto the same point
    modNew = interp1(modTime, modUwind, obsTime(iKeep));
    
    %  Error=Error1(T,P,Comment);
    tmp =  Error1(obsNew, modNew,'');
    stdUwind(iObs,1) = tmp.StdP./tmp.StdT;
    rmsUwind(iObs,1) = tmp.RMSE;
    corUwind(iObs,1) = tmp.Cor;
    
    
    % Find good observation point lying the plot window
    iKeep = find(obsTime >= tStart & obsTime <= tEnd & ~isnan(obsVwind(:,end)));
    obsNew = obsVwind(iKeep,end);
    % Map the modeled result onto the same point
    modNew = interp1(modTime, modVwind, obsTime(iKeep));
    
    %  Error=Error1(T,P,Comment);
    tmp =  Error1(obsNew, modNew,'');
    stdVwind(iObs,1) = tmp.StdP./tmp.StdT;
    rmsVwind(iObs,1) = tmp.RMSE;
    corVwind(iObs,1) = tmp.Cor;
    
end
iGo = [16 18:nObs];
stdPair(iGo) = [];
rmsPair(iGo) = [];
corPair(iGo) = [];

iGo = [14:16 18 20:nObs];
stdUwind(iGo) = [];
rmsUwind(iGo) = [];
corUwind(iGo) = [];

stdVwind(iGo) = [];
rmsVwind(iGo) = [];
corVwind(iGo) = [];

iGo = [5 16 20:nObs];
stdHwave(iGo) = [];
rmsHwave(iGo) = [];
corHwave(iGo) = [];

if 0
    stdPair = [1; stdPair(:)];
    rmsPair = [0; rmsPair(:)];
    corPair = [1; corPair(:)];
end
nStn4Pair    = length(rmsPair);

%axes(ha(2));
%subplot(1,3,2);
%[hp,ht,axl] = taylordiag_USeast(stdPair,rmsPair,corPair,'tickRMS',[0.5:0.5:1.5]);

if 0
    stdHwave = [1; stdHwave(:)];
    rmsHwave = [0; rmsHwave(:)];
    corHwave = [1; corHwave(:)];
end
nStn4Hwave   = length(rmsHwave);

%axes(ha(3));
%subplot(1,3,3);
%[hp,ht,axl] = taylordiag_USeast(stdHwave,rmsHwave,corHwave,'tickRMS',[0.5:0.5:1.5]);
stdUSeast = [1; stdTide(:); stdPair(:); stdHwave(:)];
rmsUSeast = [1; rmsTide(:); rmsPair(:); rmsHwave(:)];
corUSeast = [1; corTide(:); corPair(:); corHwave(:)];

[hp,ht,axl] = taylordiag_USeast(stdUSeast,rmsUSeast,corUSeast,'tickRMS',[0.5:0.5:1.5], ...
    'WidthCor',2, 'WidthRMS',2, 'WidthSTD', 2);


%export_fig -a1 -m1 Fig07_Taylordiagram.png
