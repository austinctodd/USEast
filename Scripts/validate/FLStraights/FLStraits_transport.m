%======================== FLStraits_transport.m ==========================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program is designed to read in results from a the US East Coast =%
%=  water age ROMS simulation and calculate the transport through the FL =%
%=  Straits at 27N, and compare to observed transport values. This       =%
%=  program requires use of the ROMS grid stretching function set_depth  =%
%=                                                                       =%
%=========================================================================%

%--- Add various libraries and paths ---%
addpath(genpath('/home/actodd/MYMATLAB/ROMS-matlab/'));

%=========================================================================%
%=                       Load ROMS data for 2004                         =%
%=========================================================================%
pdir='/he_data/he/zhigang/Project/USeast/Out_exp';
ncid=netcdf.open([pdir,'01/his_0001.nc'],'nowrite');
  Vtransform =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vtransform' ));
  Vstretching=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vstretching'));
  theta_s    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_s'    ));
  theta_b    =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_b'    ));
  hc         =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hc'         ));
  lon =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'),[195 208],[12 1]);
  lat =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'),[195 208],[12 1]);
  mask=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'));
  h   =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'),[195 208],[11 1]);
netcdf.close(ncid);

%=========================================================================%
%=                       Load ROMS data for 2004                         =%
%=========================================================================%
disp('Reading in ROMS output');

for ex=11:14
  fprintf(1,['Experiment ',sprintf('%02i',ex)]);
  if (ex ~=13)
  for i=1:6
    fprintf(1,['...Month ',sprintf('%01i',i)]);
    ncid=netcdf.open([pdir,sprintf('%02i',ex),'/his_',sprintf('%04i',i),...
                      '.nc'],'nowrite');
      zeta=double(squeeze(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'),...
                                        [195 208 0],[11 1 31])));
      v   =double(squeeze(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),...
                                        [195 207 0 0],[11 1 36 31])));
    netcdf.close(ncid);

    dist=sw_dist(lat,lon,'km').*1000;

    for j=1:size(zeta,2)
      zw=squeeze(set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,...
                      36,5,h,zeta(:,j),0));
      zw=zw(:,2:end)-zw(:,1:end-1); 
      transport(ex-10,:,j+31*(i-1))=nansum(zw.*squeeze(v(:,:,j)),2).*dist;
    end
  end
  end
  fprintf(1,'\n');
end

%=========================================================================%
%=                            Load ROMS data                             =%
%=========================================================================%
disp('Reading in FL Current data');
[FC.yr,FC.mn,FC.dy,FC.transport,FC.flag]=textread(['data/',...
    'FC_cable_transport_2013.dat'],'%4d %2d %2d %f %1d','headerlines',22);

%=========================================================================%
%=                         Calculate some stats                          =%
%=========================================================================%
for i=1:4
  if (i~=3)
  x=squeeze(sum(transport(i,:,:),2)).*1e-6;
  y=FC.transport(1:186);
  stats.rmse(i)=sqrt(mean((x-y).^2));
  stats.avg( i)=mean(x);
  p = polyfit(y,x,1);
  yfit=polyval(p,x);
  r=cov(x,y)./(std(x)*std(y));
  stats.r(i)=r(1,1); 
  stats.reg(i)=p(1);
  stats.std(i)=std(x);
  end
end
stats.FCmean=mean(y);
stats.FCstd =std(y);

%=========================================================================%
%=                               Plot data                               =%
%=========================================================================%
figure(1); clf; hold on

%--- Plot transport time series ---%
FCT=line(1:186,FC.transport(1:186));
ex1=line(1:186,squeeze(sum(transport(1,:,:),2)).*1e-6);
ex2=line(1:186,squeeze(sum(transport(2,:,:),2)).*1e-6);
%ex3=line(1:186,squeeze(sum(transport(3,:,:),2)).*1e-6);
ex4=line(1:186,squeeze(sum(transport(4,:,:),2)).*1e-6);
%ex5=line(1:186,squeeze(sum(transport(5,:,:),2)).*1e-6);

%--- Plot legend line ---%
ex1_lab=line([  3  27],[17.5 17.5]);
ex2_lab=line([ 32  56],[17.5 17.5]);
%ex3_lab=line([ 62  93],[17.5 17.5]);
ex4_lab=line([100 125],[17.5 17.5]);
%ex5_lab=line([132 157],[17.5 17.5]);
fct_lab=line([164 183],[17.5 17.5]);

%--- Write legend name ---%
ex1_nam=text(  3,16.5,['20,100 (',sprintf('%4.2f',...
                              mean(sum(transport(1,:,:),2)).*1e-6),')']);
ex2_nam=text( 32,16.5,['40,200 (',sprintf('%4.2f',...
                            mean(sum(transport(2,:,:),2)).*1e-6),')']);
%ex3_nam=text( 62,16.5,['MPDATA (',sprintf('%4.2f',...
%                            mean(sum(transport(3,:,:),2)).*1e-6),')']);
ex4_nam=text(100,16.5,['40,500 (',sprintf('%4.2f',...
                            mean(sum(transport(4,:,:),2)).*1e-6),')']);
%ex5_nam=text(131,16.5,['C2C2 (',sprintf('%4.2f',...
%                            mean(sum(transport(5,:,:),2)).*1e-6),')']);
fct_nam=text(161,16.5,[' Obs (',sprintf('%4.2f',...
                            mean(FC.transport(1:186))),')']);

%--- Set line properties ---%
set([ex1,ex1_lab],'color',[1.0 0.0 0.0],'linewidth',1.0,'linestyle','-' );
set([ex2,ex2_lab],'color',[0.0 0.7 0.0],'linewidth',1.0,'linestyle','-');
%set([ex3,ex3_lab],'color',[0.0 0.0 1.0],'linewidth',1.0,'linestyle','-');
set([ex4,ex4_lab],'color',[1.0 0.5 0.0],'linewidth',1.0,'linestyle','-' );
%set([ex5,ex5_lab],'color',[0.6 0.5 0.4],'linewidth',1.0,'linestyle','-');
%set([ex5,ex5_lab],'color',[0.0 0.0 0.0],'linewidth',1.0,'linestyle','-');
set([FCT,fct_lab],'color',[0.5 0.5 0.5],'linewidth',2.0,'linestyle','-' );

%--- Add Legends and Labels ---%
hTitle  = title ('Florida Current Transport, 2013');
hYLabel = ylabel('Transport (Sv)'                 );
hXLabel = xlabel('Days since 01 Jan 2013'         );

%--- Set Figure Properties ---%
axis([0 187 15 40]); pbaspect([2 1 1]);
set(gca,'FontName','Helvetica','Box','on','TickDir','out','Ticklength',[.01 .01],...
    'YGrid','on','XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'YMinorTick','on',...
    'XTick',0:20:180,'Linewidth',1);
set([hTitle,hXLabel,hYLabel],'FontName','Helvetica');
set([gca],'FontSize',8)
set([hXLabel,hYLabel],'FontSize',10);
set(gcf,'color','w');

%--- Adjust Title Location and Appearance ---%
hh=get(hTitle,'position');
set(hTitle,'FontSize',12,'Fontweight','bold','position',[hh(1) hh(2)+.3 hh(3)]);


return;
%=========================================================================%
%=                               Plot data                               =%
%=========================================================================%
figure(2); clf; hold on

ex1=line([1 186],[mean(sum(transport(1,:,:),2)) mean(sum(transport(1,:,:),2))].*1e-6);
ex2=line([1 186],[mean(sum(transport(2,:,:),2)) mean(sum(transport(2,:,:),2))].*1e-6);
ex3=line([1 186],[mean(sum(transport(3,:,:),2)) mean(sum(transport(3,:,:),2))].*1e-6);
ex4=line([1 186],[mean(sum(transport(4,:,:),2)) mean(sum(transport(4,:,:),2))].*1e-6);
ex5=line([1 186],[mean(sum(transport(5,:,:),2)) mean(sum(transport(5,:,:),2))].*1e-6);
FCT=line([1 186],[mean(FC.transport(1:186)) mean(FC.transport(1:186))]);

set(ex1,'color',[1.0 0.0 0.0],'linewidth',1.0,'linestyle','-' );
set(ex2,'color',[1.0 0.0 0.0],'linewidth',1.0,'linestyle','--');
set(ex3,'color',[0.0 0.0 0.0],'linewidth',1.5,'linestyle','--');
set(ex4,'color',[0.0 0.0 1.0],'linewidth',1.0,'linestyle','-' );
set(ex5,'color',[0.0 0.0 1.0],'linewidth',1.0,'linestyle','--');
set(FCT,'color',[0.0 0.0 0.0],'linewidth',1.5,'linestyle','-' );

axis([0 182 20 40]);
