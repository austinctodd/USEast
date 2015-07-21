%============================= plot_eke.m ================================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2015)                                =%
%=                                                                       =%
%=  This program is designed to read in the US East Coast ROMS grid and  =%
%=  plot starting locations of grid cells that are between 50 and 200 nm =%
%=  offshore of the Southeast.                                           =%
%=                                                                       =%
%=========================================================================%

%--- Add various libraries and paths ---%
addpath(genpath('/Users/actodd/MYMATLAB/ROMS-matlab/'));
addpath('/Users/actodd/MYMATLAB/');
addpath(genpath('/Users/actodd/MYMATLAB/m_map'));
format long g; format compact;
set(0,'defaultfigurerenderer','zbuffer');

%=========================================================================%
%=                  Set file input/output directories                    =%
%=========================================================================%
plotdir='/Users/actodd/PLOTS/USeast/temp/frames/';
datadir='/Volumes/Black_box/Data/USeast/Data/grd/';

%=========================================================================%
%=                      Load ROMS grid information                       =%
%=========================================================================%
disp('Reading in ROMS grid information');

fname=[datadir,'USeast-grid.nc'];

ncid=netcdf.open(fname,'nowrite');
 lon        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon_rho'    ),...
                                                [0 0],[402 482]);
 lat        =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat_rho'    ),...
                                                [0 0],[402 482]);
 mask       =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'mask_rho'   ),...
                                                [0 0],[402 482]);
 h          =netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'          ),...
                                                [0 0],[402 482]);
netcdf.close(ncid);

lon=lon'; lat=lat'; mask=mask'; h=h';


[sizelt,sizeln]=size(mask);
nanmask=mask./mask;

%=========================================================================%
%=                     Set different axes properties                     =%
%=========================================================================%
lats=lat(find(lat(:,1)>=30.25 & lat(:,1)<=38.5),1);
latinds=[find(lat(:,1)==lats(1)),find(lat(:,1)==lats(end))];

%=========================================================================%
% Define diagonal lines from DE/J and GA/FL borders to find points within
%=========================================================================%
mDE=(37.25-38.5)/(-71+74);
bDE=(mDE*71)+37.25;

mGA=(30-30.7)/(-77+81);
bGA=(mGA*77)+30;

%=========================================================================%
% Find lons of ocean points in that range (limit lon search to 65W)
%=========================================================================%
ln=[];
lt=[];
md=[];
for i =latinds(1):latinds(2)
    
    % Find point closest to shore
    ind=max(find(mask(i,1:308)<1));
    
    % Move Eastward
    for j=ind+1:308
        
        % find if point is within 50nm of land
        mindist=201;
        
        % Make sure it is within the north and south border lines
        if (lat(i,j)>=(mGA*lon(i,j)+bGA) && lat(i,j)<=(mDE*lon(i,j)+bDE))
            for ii=i-40:min([482,i+40])
                jc=1;
                for jj=j-40:min([402,j+40])
                    if (mask(ii,jj)<1)
                        % Find the distance of the land
                        x=[lon(i,j),lon(ii,jj)];
                        y=[lat(i,j),lat(ii,jj)];
                        [a,b]=sw_dist(y,x,'nm');

                        mindist=min([mindist,a]);
                    end
                end
            end
            if (mindist>=30 && mindist <=200)
                ln(end+1)=lon(i,j);
                lt(end+1)=lat(i,j);
                md(end+1)=mindist;
            end
        end
    end
end

return;

%=========================================================================%
%=                     Set different axes properties                     =%
%=========================================================================%
ax1pos(1,:)=[0.035 0.53 0.35 0.465];
ax1pos(2,:)=[0.340 0.53 0.35 0.465];
ax1pos(3,:)=[0.035 0.03 0.35 0.465];
ax1pos(4,:)=[0.340 0.03 0.35 0.465];
ax1pos(5,:)=[0.645 0.53 0.35 0.465];

exnames={'era';'era_hycom';'cfsr';'cfsr_hycom'};
extitle={'ERA';'ERA/HYCOM';'CFSR';'CFSR/HYCOM'};

figure(1); clf;
return;
for ex=1:4
    
  eval(['temp=squeeze(nanmean(',exnames{ex},'.dQdSST(:,:,10:12),3));']);
  
  %=======================================================================%
  %=                   Main plot panel for each exp                      =%
  %=======================================================================%
  eval(['ax',sprintf('%1i',ex),'=axes(''position'',ax1pos(',...
             sprintf('%1i',ex),',:));']);
  m_proj('lambert','long',[min(min(lon)) max(max(lon))],...
                    'lat',[min(min(lat)) max(max(lat))]);
  if (ex==1)
    m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
           'xticklabels',char({'96^{o}W';'88^{o}W';'80^{o}W';'72^{o}W';''}),...
	       'yticklabels',char({'16^{o}N ';'24^{o}N ';'32^{o}N ';'40^{o}N ';}),...
	       'fontsize',10)
    elseif (ex==2)
      m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
             'xticklabels',char({'96^{o}W';'88^{o}W';'80^{o}W';'72^{o}W';''}),...
	     'yticklabels',char({'';'';'32^{o}N ';'40^{o}N '}),...
	     'fontsize',10)
    elseif (ex==3)
      m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
             'xticklabels',char({'96^{o}W';'88^{o}W';'80^{o}W';'72^{o}W';''}),...
	     'yticklabels',char({'16^{o}N ';'24^{o}N ';'32^{o}N ';'40^{o}N ';}),...
	     'fontsize',10)
    else
      m_grid('tickdir','on','xtick',[-96:8:-64],'ytick',[16:8:40],...
             'xticklabels',char({'96^{o}W';'88^{o}W';'80^{o}W';'72^{o}W';'64^{o}W'}),...
	     'yticklabels',char({'';'';'32^{o}N ';'40^{o}N ';}),...
	     'fontsize',10)
    end
    hold on;
    m_contour(lon,lat,mask,[0 0],'k','linewidth',1.5);
    pcc=m_pcolor(lon,lat,-temp.*mask./mask); shading flat; 
    caxis([8 30]); 

    m_text(-96.5,43,extitle{ex},'fontsize',10,'fontweight','bold');
end

%=========================================================================%
%=                             Add Colorbar                              =%
%=========================================================================%
cax=axes('position',[0.7 0.05 0.1 0.65]); axis off
caxis([8.0 30])
hCbar     =colorbar('location','East','position',[0.7 0.05 0.1 0.67],...
                     'yaxislocation','right');
hCbarT=text(1.65,0.75,'dQ/dSST (Wm^{-2}K^{-1})','fontsize',12,...
                      'fontname','helvetica','fontweight','bold',...
                      'rotation',270);
set(hCbar,'linewidth',1,'xtick',[0:5:30],'ticklength',[.02 .02],...
          'fontsize',10,'fontname','helvetica','fontweight','bold');


%=========================================================================%
%=              Make set of text at top for model run info               =%
%=========================================================================%
tax=axes('position',[0.7 0.27 0.275 0.2],'xtick',[],'ytick',[]); box on;
set(tax,'position',[0.67 0.75 0.225 0.225]);
set(tax,'linewidth',1); 
hmodelname=text(0.03,0.9,'US East Coast ROMS','fontweight','bold',...
    'fontsize',12);
hNCSU     =text(0.065,0.75,'NC State Univeristy (2014)','fontsize',10);             
hmodelres =text(0.45,0.60,'\Delta x,\Delta y = 1/12^{o}','fontsize',10);          
hmodelvert=text(0.76,0.46,'k = 36','fontsize',10);
hmodeldate=text(0.975,0.2,'Oct-Dec dQ/dSST',...
                  'fontsize',14,'fontweight','bold',...
                  'horizontalalignment','right');

set(gcf,'color','w');
%axes(ax5); 
axes(ax2); axes(ax3); axes(ax1);
wysiwyg