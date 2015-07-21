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
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));

%=========================================================================%
%=                      Set list of all station IDs                      =%
%=========================================================================%
load('/he_data/he/actodd/DATA/tides/GOM/GOM_harmonics.mat');
tide.M2amp=M2amp; tide.M2phase=M2phase;
tide.S2amp=S2amp; tide.S2phase=S2phase;
tide.O1amp=O1amp; tide.K1phase=O1phase;
tide.K1amp=K1amp; tide.O1phase=K1phase;

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
stnName={'Loggerhead Key';'Vaca Key';'Key West';'Fort Myers';...
         'Port Manatee';'Clearwater';'Johns Island';'Shell Island';...
         'St. Marks Lighthouse';'Turkey Point';'Apalachicola';...
         'Panama City';'Navarre Beach';'Pensacola';...
         'Dauphin Island Hydro';'Horn Island, MS';'Gulfport Harbor';...
         'Waveland';'South Pass';'Grand Isle';'Port Fourchon';...
         'Eugene Island';'Freswater Canal Locks';'Calcasieu';...
         'Sabine offshore';'RolloverPass';'Galveston';'USCG Freeport';...
         'Port O''Connor';'Port Aransas, Caldwell Pier';...
         'Packery Channel';'Corpus Christi';'Padre Island';...
         'Madero, Tampico';'Naples';'St. Petersburg';'Tuckers Island';...
         'Mangrove Point';'Ozello North';'Cedar Key';'Shell Point';...
         'Weeks Bay';'Pascagoula Point';'Coast Guard sector Mobile';...
         'Biloxi';'Bay Waveland Yacht Club';'Shell Beach';...
         'South West Pass';'Texaco Dock';'Lawma, Amerada Pass';...
         'Cypremort Point';'Galveston Bay entrance';...
         'Port Ingleside, Corpus Christi Bay';'Port Isabel'};
      
      
stnLon=-(stnLonD+stnLonM./60);
stnLat= stnLatD+stnLatM./60;

ss.ID=stnID; ss.lat=stnLat; ss.lon=stnLon;
      
nSta=length(stnID);

%=========================================================================%
%=                   Set ROMS file name and load data                    =%
%=========================================================================%
roms.f='/gpfs_share/rhe/actodd/USeast-tide/output/LTP/orig/useast_his.nc';
roms.lon =nc_varget(roms.f,'lon_rho' ,[80 9],[183 177]);
roms.lat =nc_varget(roms.f,'lat_rho' ,[80 9],[183 177]);
roms.mask=nc_varget(roms.f,'mask_rho',[80 9],[183 177]);

data_dir='/gpfs_share/rhe/actodd/USeast-tide/output/';

%=========================================================================%
%=                 Load LTP/OBC harmonic analysis data                   =%
%=========================================================================%
load([data_dir,'LTP_OBC/newtide/GOM_harmonic_analysis.mat']);
ltpobc.M2amp=M2amp; ltpobc.M2phase=M2phase;
ltpobc.S2amp=S2amp; ltpobc.S2phase=S2phase;
ltpobc.O1amp=O1amp; ltpobc.K1phase=O1phase;
ltpobc.K1amp=K1amp; ltpobc.O1phase=K1phase;

%=========================================================================%
%=                   Load LTP harmonic analysis data                     =%
%=========================================================================%
load([data_dir,'LTP/newtide/GOM_harmonic_analysis.mat']);
ltp.M2amp=M2amp; ltp.M2phase=M2phase;
ltp.S2amp=S2amp; ltp.S2phase=S2phase;
ltp.O1amp=O1amp; ltp.K1phase=O1phase;
ltp.K1amp=K1amp; ltp.O1phase=K1phase;

%=========================================================================%
%=                   Load OBC harmonic analysis data                     =%
%=========================================================================%
load([data_dir,'OBC/newtide/GOM_harmonic_analysis.mat']);
obc.M2amp=M2amp; obc.M2phase=M2phase;
obc.S2amp=S2amp; obc.S2phase=S2phase;
obc.O1amp=O1amp; obc.K1phase=O1phase;
obc.K1amp=K1amp; obc.O1phase=K1phase;
clear M2amp M2phase S2amp S2phase O1amp O1phase K1amp K1phase

%=========================================================================%
%=                   Set Tide file name and load data                    =%
%=========================================================================%
xint=0.0125; xdist=(1-(xint*20))/4;
yint=0.0075; ydist=(1-(yint*22))/18;

%=========================================================================%
%=                   Set Tide file name and load data                    =%
%=========================================================================%
count=1; figure(1); clf; figure(2); clf;
for i=1:nSta
  
  if tide.M2amp(i)>0
    [xind,yind]=find_closest_roms_point(roms,ss,stnID(i));
    
    figure(1);
    y1=yint+yint*(count-1)+ydist*(count-1);
    eval(['fig1axes'  ,sprintf('%02i',count),'_title',...
             '=axes(''position'',[',sprintf('%8.6f',xint),',',...
                                    sprintf('%8.6f',y1   ),',',...
                                    sprintf('%8.6f',xdist),',',...
                                    sprintf('%8.6f',ydist),']);']);
    axis off;
    text(0.45,.5,stnName{i},'fontsize',9,'horizontalalignment','center');
       
    figure(2);
    eval(['fig2axes'  ,sprintf('%02i',count),'_title',...
             '=axes(''position'',[',sprintf('%8.6f',xint),',',...
                                    sprintf('%8.6f',y1   ),',',...
                                    sprintf('%8.6f',xdist),',',...
                                    sprintf('%8.6f',ydist),']);']);
    axis off;
    text(0.45,.5,stnName{i},'fontsize',9,'horizontalalignment','center');
    
    for const=1:4
      
      %--- Make axes for plotting tide amplitude ---%
      figure(1);
      x1=xint*15+xint*(const-1)+xdist*(const-1);
      y1=yint+yint*(count-1)+ydist*(count-1);
      eval(['fig1axes'  ,sprintf('%02i',count),'_',sprintf('%1i',const),...
               '=axes(''position'',[',sprintf('%8.6f',x1   ),',',...
                                      sprintf('%8.6f',y1   ),',',...
                                      sprintf('%8.6f',xdist),',',...
                                      sprintf('%8.6f',ydist),']);']);
 
      %--- Plot tidal amp like Gouillon et al. 2010 (Fig. 7) ---%
      if (const==1)
        hold on;
        bar(1,  tide.M2amp(i        ),'Facecolor',[0.2 0.2 0.2]);
        bar(2,ltpobc.M2amp(yind,xind),'Facecolor',[0.4 0.4 0.4]); 
        bar(3,   obc.M2amp(yind,xind),'Facecolor',[0.6 0.6 0.6]);
        bar(4,   ltp.M2amp(yind,xind),'Facecolor',[0.8 0.8 0.8]);
        axis([0.25 4.75 0 0.40]);
      elseif (const==2)
        hold on;
        bar(1,  tide.S2amp(i        ),'Facecolor',[0.2 0.2 0.2]);
        bar(2,ltpobc.S2amp(yind,xind),'Facecolor',[0.4 0.4 0.4]);
        bar(3,   obc.S2amp(yind,xind),'Facecolor',[0.6 0.6 0.6]);
        bar(4,   ltp.S2amp(yind,xind),'Facecolor',[0.8 0.8 0.8]);
        axis([0.25 4.75 0 0.18]);
      elseif (const==3)
        hold on;
        bar(1,  tide.O1amp(i        ),'Facecolor',[0.2 0.2 0.2]);
        bar(2,ltpobc.O1amp(yind,xind),'Facecolor',[0.4 0.4 0.4]);
        bar(3,   obc.O1amp(yind,xind),'Facecolor',[0.6 0.6 0.6]);
        bar(4,   ltp.O1amp(yind,xind),'Facecolor',[0.8 0.8 0.8]);
        axis([0.25 4.75 0 0.22]); 
      else
        hold on;
        bar(1,  tide.K1amp(i        ),'Facecolor',[0.2 0.2 0.2]);
        bar(2,ltpobc.K1amp(yind,xind),'Facecolor',[0.4 0.4 0.4]);
        bar(3,   obc.K1amp(yind,xind),'Facecolor',[0.6 0.6 0.6]);
        bar(4,   ltp.K1amp(yind,xind),'Facecolor',[0.8 0.8 0.8]);
        axis([0.25 4.75 0 0.24]);
      end
      set(gca,'box','on','xtick',[],'yticklabel',[]);
   
      %--- Make axes for plotting tide phase ---%
      figure(2);
      x1=xint*15+xint*(const-1)+xdist*(const-1);
      y1=yint+yint*(count-1)+ydist*(count-1);
      eval(['fig2axes'  ,sprintf('%02i',count),'_',sprintf('%1i',const),...
               '=axes(''position'',[',sprintf('%8.6f',x1   ),',',...
                                      sprintf('%8.6f',y1   ),',',...
                                      sprintf('%8.6f',xdist),',',...
                                      sprintf('%8.6f',ydist),']);']);
      plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k','linewidth',1);hold on;
                                  
      %--- Plot tidal phase like Gouillon et al. 2010 (Fig. 7) ---%
      if (const==1)
        plot([0 cosd(  tide.M2phase(i        ))],...
             [0 sind(  tide.M2phase(i        ))],'k'  ,'linewidth',0.5);
        plot([0 cosd(ltpobc.M2phase(yind,xind))],...
             [0 sind(ltpobc.M2phase(yind,xind))],'k:','linewidth',0.5);
        plot([0 cosd(   obc.M2phase(yind,xind))],...
             [0 sind(   obc.M2phase(yind,xind))],'k:','linewidth',0.5,...
             'color',[.5 .5 .5]);
        plot([0 cosd(   ltp.M2phase(yind,xind))],...
             [0 sind(   ltp.M2phase(yind,xind))],'k'  ,'linewidth',0.5,...
             'color',[.5 .5 .5]);
        axis equal; axis([-1.05 1.05 -1.05 1.05]); axis off;
      elseif (const==2)
        plot([0 cosd(  tide.S2phase(i        ))],...
             [0 sind(  tide.S2phase(i        ))],'k'  ,'linewidth',0.5);
        plot([0 cosd(ltpobc.S2phase(yind,xind))],...
             [0 sind(ltpobc.S2phase(yind,xind))],'k:','linewidth',0.5);
        plot([0 cosd(   obc.S2phase(yind,xind))],...
             [0 sind(   obc.S2phase(yind,xind))],'k:','linewidth',0.5,...
             'color',[.5 .5 .5]);
        plot([0 cosd(   ltp.S2phase(yind,xind))],...
             [0 sind(   ltp.S2phase(yind,xind))],'k'  ,'linewidth',0.5,...
             'color',[.5 .5 .5]);
        axis equal; axis([-1.05 1.05 -1.05 1.05]); axis off;
      elseif (const==3)
        plot([0 cosd(  tide.O1phase(i        ))],...
             [0 sind(  tide.O1phase(i        ))],'k'  ,'linewidth',0.5);
        plot([0 cosd(ltpobc.O1phase(yind,xind))],...
             [0 sind(ltpobc.O1phase(yind,xind))],'k:','linewidth',0.5);
        plot([0 cosd(   obc.O1phase(yind,xind))],...
             [0 sind(   obc.O1phase(yind,xind))],'k:','linewidth',0.5,...
             'color',[.5 .5 .5]);
        plot([0 cosd(   ltp.O1phase(yind,xind))],...
             [0 sind(   ltp.O1phase(yind,xind))],'k'  ,'linewidth',0.5,...
             'color',[.5 .5 .5]);
        axis equal; axis([-1.05 1.05 -1.05 1.05]); axis off;
      else
        plot([0 cosd(  tide.K1phase(i        ))],...
             [0 sind(  tide.K1phase(i        ))],'k'  ,'linewidth',0.5);
        plot([0 cosd(ltpobc.K1phase(yind,xind))],...
             [0 sind(ltpobc.K1phase(yind,xind))],'k:','linewidth',0.5);
        plot([0 cosd(   obc.K1phase(yind,xind))],...
             [0 sind(   obc.K1phase(yind,xind))],'k:','linewidth',0.5,...
             'color',[.5 .5 .5]);
        plot([0 cosd(   ltp.K1phase(yind,xind))],...
             [0 sind(   ltp.K1phase(yind,xind))],'k'  ,'linewidth',0.5,...
             'color',[.5 .5 .5]);
        axis equal; axis([-1.05 1.05 -1.05 1.05]); axis off;
      end
      set(gca,'box','off','xtick',[],'yticklabel',[]);
    end
    count=count+1;
  end
end

figure(1);
eval(['axes(fig1axes',sprintf('%02i',count-1),'_1',');']);
titM2=title('M2','fontsize',10','fontweight','bold');
pos=get(titM2,'position');
set(titM2,'position',[pos(1) pos(2)*0.85 pos(3)]);
clear pos;

eval(['axes(fig1axes',sprintf('%02i',count-1),'_2',');']);
titS2=title('S2','fontsize',10','fontweight','bold');
pos=get(titS2,'position');
set(titS2,'position',[pos(1) pos(2)*0.85 pos(3)]);
clear pos;

eval(['axes(fig1axes',sprintf('%02i',count-1),'_3',');']);
titO1=title('O1','fontsize',10','fontweight','bold');
pos=get(titO1,'position');
set(titO1,'position',[pos(1) pos(2)*0.85 pos(3)]);
clear pos;

eval(['axes(fig1axes',sprintf('%02i',count-1),'_4',');']);
titK1=title('K1','fontsize',10','fontweight','bold');
pos=get(titK1,'position');
set(titK1,'position',[pos(1) pos(2)*0.85 pos(3)]);

figure(2);
eval(['axes(fig2axes',sprintf('%02i',count-1),'_1',');']);
titM2=title('M2','fontsize',10','fontweight','bold');
pos=get(titM2,'position');
set(titM2,'position',[pos(1) pos(2)*0.85 pos(3)]);
clear pos;

eval(['axes(fig2axes',sprintf('%02i',count-1),'_2',');']);
titS2=title('S2','fontsize',10','fontweight','bold');
pos=get(titS2,'position');
set(titS2,'position',[pos(1) pos(2)*0.85 pos(3)]);
clear pos;

eval(['axes(fig2axes',sprintf('%02i',count-1),'_3',');']);
titO1=title('O1','fontsize',10','fontweight','bold');
pos=get(titO1,'position');
set(titO1,'position',[pos(1) pos(2)*0.85 pos(3)]);
clear pos;

eval(['axes(fig2axes',sprintf('%02i',count-1),'_4',');']);
titK1=title('K1','fontsize',10','fontweight','bold');
pos=get(titK1,'position');
set(titK1,'position',[pos(1) pos(2)*0.85 pos(3)]);



