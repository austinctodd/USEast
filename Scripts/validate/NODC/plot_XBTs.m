%============================= plot_XBTs.m ===============================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program compare results from a the US East Coast water age ROMS =%
%=  simulation to a whole set of XBT casts from the NODC datasets.       =%
%=                                                                       =%
%=========================================================================%

addpath('/home/actodd/MYMATLAB/');
%load /home/actodd/ROMS-utils/USeast-age/validate/XBT_ROMS_data.mat
%cruisenames=unique(XBT.cruise);

for i=1:size(cruisenames,1)

  if length(XBT.cruisetemp{i})>=10

    figure(1); clf;
    ax1  =axes('position',[0.05 0.520 0.80 0.220]);
    ax2  =axes('position',[0.05 0.295 0.80 0.220]);
    cbAx1=axes('position',[0.86 0.295 0.04 0.445]); axis off;
  
    %--- Initialize axis limits ---%
    zbot=1e36; ztop=-1e36;
    smin=1e36; smax=-1e36;
    tmin=1e36; tmax=-1e36;
    xmin=1e36; xmax=-1e36;
  
    for j=1:length(XBT.cruisetemp{i})
      %===================================================================%
      %=                 Observational temperature plot                  =%
      %===================================================================%
      axes(ax1);hold on;
           
      %-- Plot profile --%
      scatter( XBT.cruisetime(i,j).*ones(length( XBT.cruisetemp{i}{j}),1),...
              -XBT.cruisez{i}{j}{:},20, t,'filled');
  
      %===================================================================%
      %=                     Model temperature plot                      =%
      %===================================================================%
      axes(ax2);hold on;
      
      %-- Plot profile --%
      scatter( XBT.cruisetime(i,j).*ones(length(roms.cruisetemp{i}{j}),1),...
              roms.cruisez{i}{j},20,roms.cruisetemp{i}{j},'filled');
    
      %===================================================================%
      %=                        Find axis limits                         =%
      %===================================================================%
      zbot=min([zbot nanmin(-XBT.cruisez{i}{j}{:})]);% min(roms.cruisez{i}{j})]);
      ztop=max([ztop nanmax(-XBT.cruisez{i}{j}{:})]);% max(roms.cruisez{i}{j})]);
      tmin=min([tmin nanmin(XBT.cruisetemp{i}{j}) nanmin(roms.cruisetemp{i}{j})]);
      tmax=max([tmax nanmax(XBT.cruisetemp{i}{j}) nanmax(roms.cruisetemp{i}{j})]);
    
    end
    %=====================================================================%
    %=               Adjust observational temperature plot               =%
    %=====================================================================%
    axes(ax1);
    caxis([tmin tmax]);
    datetick;
    xlims=get(gca,'xlim');
    xlims=[min([XBT.cruisetime(i,1:length(XBT.cruisetemp{i})) ...
               roms.cruisetime(i,1:length(XBT.cruisetemp{i}))]) ...
           max([XBT.cruisetime(i,1:length(XBT.cruisetemp{i})) ...
               roms.cruisetime(i,1:length(XBT.cruisetemp{i}))])];
    xdist=xlims(2)-xlims(1);
    axis([xlims(1)-xdist/30 xlims(2)+xdist/30 zbot-5 ztop+1]);
    box on; 
    set(ax1,'linewidth',1,'fontsize',8,'xticklabel','');
    TobsText=text(xlims(1),zbot-zbot/14,['XBT : ',cruisenames(i,1:10)],...
                  'HorizontalAlignment','left',...
		  'fontsize',10,'fontweight','bold');

    %=====================================================================%
    %=                   Adjust model temperature plot                   =%
    %=====================================================================%
    axes(ax2);
    caxis([tmin tmax]);
    datetick;
    axis([xlims(1)-xdist/30 xlims(2)+xdist/30 zbot-5 ztop+1]);
    set(ax2,'linewidth',1,'fontsize',8);
    box on; 
    TmodText=text(xlims(1),zbot-zbot/14,'US East ROMS',...
                  'HorizontalAlignment','left',...
		          'fontsize',10,'fontweight','bold');
		  
    %=====================================================================%
    %=                     Make temperature colorbar                     =%
    %=====================================================================%
    axes(cbAx1);
    cb1=colorbar('position',[0.86 0.295 0.04 0.445]);
    caxis([tmin tmax]);
    set(cb1,'linewidth',1,'fontsize',8);
    cbText1=text(2,.5,'Temperature ( ^{o}C)','fontsize',12,'rotation',270,...
                    'HorizontalAlignment','center');
                
    %=====================================================================%
    %=                     Export figure to file                         =%
    %=====================================================================%
    set(gcf,'color','w');
    wysiwyg;
    eval(['export_fig /home/actodd/PLOTS/USeast/validate/XBT/',...
           cruisenames(i,1:8),'.png -png -r150 -painters']);
    
    %=====================================================================%
    %=          Make figure with cruise track and model points           =%
    %=====================================================================%
    figure(2); clf;
    
    %--- Contour land mask from model and model depths ---%
    contour(roms.lon,roms.lat,roms.mask,[0 0],'k','linewidth',2); hold on
    contour(roms.lon,roms.lat,roms.h,[50 100 250 500 7500 1000:1000:4000],...
            'color',[.4 .4 .4],'linewidth',0.75);

    %--- Plot ROMS points as plus signs ---%
    plot(roms.cruiselon( i,1:length(XBT.cruisetemp{i})),...
         roms.cruiselat( i,1:length(XBT.cruisetemp{i})),'m+');
   
    %--- Plot cruise track and XBT cast sites ---%
    plot(XBT.cruiselon( i,1:length(XBT.cruisetemp{i})),...
         XBT.cruiselat( i,1:length(XBT.cruisetemp{i})),'k','linewidth',1);
    scatter(XBT.cruiselon( i,1:length(XBT.cruisetemp{i})),...
            XBT.cruiselat( i,1:length(XBT.cruisetemp{i})),25,...
            XBT.cruisetime(i,1:length(XBT.cruisetemp{i})),'filled');
    caxis([XBT.cruisetime(i,1) XBT.cruisetime(i,length(XBT.cruisetemp{i}))])  
    
    %--- Add colorbar for time scale --%
    hCbar=colorbar;
    if (diff(XBT.cruisetime(i,[1 length(XBT.cruisetemp{i})]))>2)
      cbartickstyle=6;
    else
      cbartickstyle=15;    
    end
    set(hCbar,'linewidth',1,'fontsize',8,'yticklabel',...
        datestr(get(hCbar,'ytick'),cbartickstyle));
    
    %--- Set axis and figure properties ---%
    axis xy; axis image;
    axis([min(XBT.cruiselon(i,1:length(XBT.cruisetemp{i})))-1 ...
          max(XBT.cruiselon(i,1:length(XBT.cruisetemp{i})))+1 ...
          min(XBT.cruiselat(i,1:length(XBT.cruisetemp{i})))-1 ...
          max(XBT.cruiselat(i,1:length(XBT.cruisetemp{i})))+1])
    
    set(gcf,'color','w');
    set(gca,'linewidth',1,'fontsize',10,'fontname','Helvetica');
    
    %--- Set title ---%
    title(['Cruise ',cruisenames(i,1:8)],'fontsize',12,'fontname',...
          'Helvetica','fontweight','bold');
      
    %--- Export figure ---%
    wysiwyg;
    eval(['export_fig /home/actodd/PLOTS/USeast/validate/XBT/',...
           cruisenames(i,1:8),'_track.png -png -r150 -painters']);

  end
end
    
