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
%=                         User-defined indexes                          =%
%=========================================================================%
tStart = datenum(2004,01,01);
tEnd   = datenum(2006,01,01);

%=========================================================================%
%=                   Set ROMS file name and load data                    =%
%=========================================================================%
roms.dir ='/gpfs_share/actodd/USeast-age/output/';
roms.f   =[roms.dir,'2004/useast_his.nc'];
ncid=netcdf.open(roms.f,'nowrite');
  roms.lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho' ))';
  roms.lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho' ))';
  roms.mask=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'))';
  roms.time=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'));
  roms.iniz=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'),...
                                               [0 0 0],[402 482 1])';
netcdf.close(ncid);
  
%roms.time=nc_varget(roms.f,'ocean_time',481,736);
%roms.date=roms.time./86400+datenum(2004,1,1);

%=========================================================================%
%=                   Set CFSR file name and load data                    =%
%=========================================================================%
cfsr.f='/he_data/he/actodd/ROMS-age/Data/frc/USeast-frc.nc';
ncid=netcdf.open(cfsr.f,'nowrite');
 cfsr.time(    1:367      )=netcdf.getVar(ncid,netcdf.inqVarID(ncid,...
                                          'time'),0,367,24);
 cfsr.time(367+1:367+365*16-1)=netcdf.getVar(ncid,netcdf.inqVarID(ncid,...
                                          'time'),24*366+3,365*16-1,3);
 cfsr.date=cfsr.time+datenum(2004,1,1);

 cfsr.lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
 cfsr.lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
netcdf.close(ncid);

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
%=                   Set Tide file name and load data                    =%
%=========================================================================%
for i=1:nSta
  disp(List4sID(i))
  tide.f=['/he_data/he/actodd/DATA/tides/',sprintf('%7i',List4sID(i)),...
          '.dat'];
  tide.data=load(tide.f);
  tide.date=datenum(tide.data(:,1),tide.data(:,2),tide.data(:,3),...
                    tide.data(:,4),tide.data(:,5),tide.data(:,5).*0);
  
  [xind,yind]=find_closest_roms_point(roms,ss,List4sID(i));
  %--- Load in 2004 data ---%
  ncid=netcdf.open([roms.dir,'2004/useast_his.nc'],'nowrite');
    roms.time(1:367)=netcdf.getVar(ncid,netcdf.inqVarID(ncid,...
                                           'ocean_time'),0,367);
    roms.zeta(1:367)=netcdf.getVar(ncid,netcdf.inqVarID(ncid,...
                                           'zeta'),[xind-1 yind-1   0],...
                                                   [     1      1 367]);
  netcdf.close(ncid);
  
  %--- Load in 2005 data ---%
  ncid=netcdf.open([roms.dir,'2005/useast_his.nc'],'nowrite');
    roms.time(367+1:367+365*16-1)=netcdf.getVar(ncid,netcdf.inqVarID(ncid,...
                                           'ocean_time'),1,365*16-1);
    roms.zeta(367+1:367+365*16-1)=netcdf.getVar(ncid,netcdf.inqVarID(ncid,...
                                           'zeta'),[xind-1 yind-1        1],...
                                                   [     1      1 365*16-1]);
  netcdf.close(ncid);
  
  roms.zeta=roms.zeta-roms.iniz(yind,xind);
  roms.date=roms.time./86400+datenum(2004,1,1);
  
  [xind,yind]=find_closest_cfsr_point(cfsr,ss,List4sID(i));
  ncid=netcdf.open(cfsr.f,'nowrite');
    cfsr.pres(    1:367      )=netcdf.getVar(ncid,netcdf.inqVarID(ncid,...
                                             'Pair'),[xind-1 yind-1 0],...
                                             [1 1 367],[1 1 24]); 
    cfsr.pres(367+1:367+365*16-1)=netcdf.getVar(ncid,netcdf.inqVarID(ncid,...
                                           'Pair'),[xind-1 yind-1 24*366+3],...
                                            [1 1 365*16-1],[1 1 3]);
  netcdf.close(ncid);   

  
  %--- Find correlation coefficients ---%
  x=roms.zeta-(9.948/1000).*(cfsr.pres./100-1013.3);
  %x=x(61:end);
  for j=1:length(roms.date);
     ind=find(tide.date==roms.date(j));
     if (length(ind)>0)
         y(j)=tide.data(ind,6);
     else
         y(j)=nan;
     end
  end
  
  ind=find(tide.date==datenum(2005,1,1));
  %y=tide.data([1:30*8:ind ind+30:30:end],6);
  a=find(isnan(y)==0); x=x(a); y=y(a);
  a=find(isnan(x)==0); x=x(a); y=y(a);
  r=cov(x,y)./(std(x)*std(y));  
 
  roms.r(  i)=r(1,2);
  roms.reg(i)=std(x)./std(y);
  roms.std(i)=std(x);
  tide.std(i)=std(y);
end

r=roms.r; reg=roms.reg; rstd=roms.std; tstd=tide.std;
save /he_data/he/actodd/DATA/tides/roms_tide_comp.mat r reg rstd tstd -MAT

  %figure(1); clf;
  %subplot(3,1,1)
  %plot(roms.date,roms.zeta-(9.948/1000).*(cfsr.pres./100-1013.3),'k','linewidth',1); hold on
  %plot(tide.date([1:30*8:ind ind+30:30:end]),tide.data([1:30*8:ind ind+30:30:end],6),'r','linewidth',1);
  %datetick
  %axis([datenum(2005,1,1),datenum(2005,2,1),...
  %       -max(abs(roms.zeta))-0.5*nanstd(roms.zeta),...
  %        max(abs(roms.zeta))+0.5*nanstd(roms.zeta)]);
  %title(['Station ',sprintf('%7i',List4sID(i)),' (R=',...
  %                  sprintf('%4.2f',r(1,2)),')'],'fontsize',12,...
  %                  'fontweight','bold');
  %set(gca,'xtick',datenum(2005,1,[1 5 10 15 20 25 32]),'xticklabel',...
  %        datestr(datenum(2005,1,[1 5 10 15 20 25 32]),6));
  %
  %subplot(3,1,2)
  %plot(roms.date,roms.zeta,'k','linewidth',1); hold on
  %plot(tide.date,tide.data(:,6),'r','linewidth',1);
  %datetick;
  %axis([datenum(2004,4,1),datenum(2004,5,1),...
  %     -max(abs(roms.zeta))-0.5*std(roms.zeta),...
  %      max(abs(roms.zeta))+0.5*std(roms.zeta)]);
  
  %subplot(3,1,3)
  %plot(roms.date,roms.zeta,'k','linewidth',1); hold on
  %plot(tide.date,tide.data(:,6),'r','linewidth',1);
  %datetick;
  %axis([datenum(2004,5,1),datenum(2004,6,1),...
  %     -max(abs(roms.zeta))-0.5*std(roms.zeta),...
  %      max(abs(roms.zeta))+0.5*std(roms.zeta)]);
 %
 % set(gcf,'color','w'); 
 %drawnow;

 %eval(['export_fig /home/actodd/PLOTS/USeast/tide_gauges/',...
 %        sprintf('%7i',List4sID(i)),'.png -png -r150 -painters']);
%end




