%============================== AVISO_EKE.m ==============================%
%=                                                                       =%
%=  Written by Austin C. Todd, NCSU (2015) for personal use              =%
%=                                                                       =%
%=  Program loads AVISO Geostrophic velocity anomolies and calculates    =%
%=  the EKE over the US East ROMS model domain.                          =%
%=                                                                       =%
%=========================================================================%

addpath(genpath('/Users/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/Users/actodd/MYMATLAB/');
addpath(genpath('/Users/actodd/MYMATLAB/m_map'));
format long g; format compact;

%=========================================================================%
%=                   Set file input/output directories                   =%
%=========================================================================%
AVISOdir ='/Volumes/Black_box/Data/AVISO/SLA/';
ROMSdir  ='/Volumes/Black_box/Data/USeast-age/output/clim/averages/';

%=========================================================================%
%=                   Set file input/output directories                   =%
%=========================================================================%
AVISOurl =['http://opendap.aviso.altimetry.fr/thredds/dodsC/',...
           'dataset-duacs-dt-global-allsat-msla-h?sla'];
latinds  ='[383:1:548]';
loninds  ='[1043:1:1204]';
fname=[AVISOurl,'[0:1:0]',latinds,loninds,',crs'];
loaddap(fname);

%--------------------------------------------------------------------------
% Declare constants used for calculations
%--------------------------------------------------------------------------
g       = -9.81;
for x=1:length(sla.lon)
  [dy(:,x),tr]=sw_dist(sla.lat,ones(length(sla.lat),1)*(sla.lon(x)-360),'km');
  f(:,x)=sw_f(sla.lat);
end
dy=dy.*1000;
for y=1:length(sla.lat)
  [dx(y,:),tr]=sw_dist(ones(length(sla.lon),1)*(sla.lat(y)),sla.lon-360,'km');
end
fx=(f(:,1:end-1)+f(:,2:end))./2.0;
fy=(f(1:end-1,:)+f(2:end,:))./2.0;
dx=dx.*1000;

%=========================================================================%
%=                           Loop through files                          =%
%=========================================================================%
EKE=zeros(165,161);
count=0;
return;
for tdim = 0:7669

  stdim=sprintf('%i',tdim);
  disp(['Time ',sprintf('%04i',tdim+1),'/7670']);
  
  %--- Construct url and load data ---$  
  fname=[AVISOurl,'[',stdim,':1:',stdim,']',latinds,loninds,',crs'];
  loaddap(fname);
     
   %--- Data management ---%
  sla.sla=double(sla.sla).*.0001;
  sla.sla(find(sla.sla<=-21474))=nan;

  ug  =-(g./fy./dy).*(sla.sla(2:end,:)-sla.sla(1:end-1,:));
  vg  = (g./fx./dx).*(sla.sla(:,2:end)-sla.sla(:,1:end-1));

  %--- Calculate EKE ---% 
  EKE = EKE + 0.5*(((ug(:,1:end-1)+ug(:,2:end))./2.0).^2 + ...
                   ((vg(1:end-1,:)+vg(2:end,:))./2.0).^2);
  count=count+1;
end


