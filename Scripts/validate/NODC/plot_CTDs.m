%============================= plot_CTDs.m ===============================%
%=                                                                       =%
%=  Written by Austin C Todd, NCSU (2014)                                =%
%=                                                                       =%
%=  This program compare results from a the US East Coast water age ROMS =%
%=  simulation to a whole set of CTD casts from the NODC datasets.       =%
%=                                                                       =%
%=========================================================================%

addpath('/home/actodd/MYMATLAB/');
load /home/actodd/ROMS-utils/USeast-age/validate/CTD_ROMS_data.mat
cruisenames=unique(CTD.cruise);

for i=1:length(unique(CTD.cruise))

  if length(CTD.cruisetemp{i})>=10

    figure(1); clf;
    ax1  =axes('position',[0.05 0.770 0.80 0.220]);
    ax2  =axes('position',[0.05 0.545 0.80 0.220]);
    ax3  =axes('position',[0.05 0.265 0.80 0.220]);
    ax4  =axes('position',[0.05 0.040 0.80 0.220]);
    cbAx1=axes('position',[0.86 0.545 0.04 0.445]); axis off;
    cbAx2=axes('position',[0.86 0.040 0.04 0.445]); axis off;
  
    %--- Initialize axis limits ---%
    zbot=1e36; ztop=-1e36;
    smin=1e36; smax=-1e36;
    tmin=1e36; tmax=-1e36;
    xmin=1e36; xmax=-1e36;
  
    for j=1:length(CTD.cruisetemp{i})
      %===================================================================%
      %=                 Observational temperature plot                  =%
      %===================================================================%
      axes(ax1);hold on;
      
      %--- See if extraneous values exist ---%
      t=CTD.cruisetemp{i}{j};
      a=find(t<-40 | t>50);
      if (length(a)>0); t(a)=nan;  end
      
      %-- Plot profile --%
      scatter( CTD.cruisetime(i,j).*ones(length( CTD.cruisetemp{i}{j}),1),...
              -CTD.cruisez{i}{j}{:},20, t,'filled');
  
      %===================================================================%
      %=                     Model temperature plot                      =%
      %===================================================================%
      axes(ax2);hold on;
      
      %-- Plot profile --%
      scatter( CTD.cruisetime(i,j).*ones(length(roms.cruisetemp{i}{j}),1),...
              roms.cruisez{i}{j},20,roms.cruisetemp{i}{j},'filled');
    
      %===================================================================%
      %=                   Observational salinity plot                   =%
      %===================================================================%
      axes(ax3);hold on;
      if (length(CTD.cruisesalt{i}{j})>0)
            
        %--- See if extraneous values exist ---%
        t=CTD.cruisesalt{i}{j};
        a=find(t<-40 | t>50);
        if (length(a)>0); t(a)=nan;  end
      
        %-- Plot profile --%
        scatter( CTD.cruisetime(i,j).*ones(length( CTD.cruisesalt{i}{j}),1),...
                -CTD.cruisez{i}{j}{:},20, t,'filled');
      end
          
      %===================================================================%
      %=                     Model temperature plot                      =%
      %===================================================================%
      axes(ax4);hold on;
       
      %-- Plot profile --%
      scatter( CTD.cruisetime(i,j).*ones(length(roms.cruisesalt{i}{j}),1),...
              roms.cruisez{i}{j},20,roms.cruisesalt{i}{j},'filled');
    
      %===================================================================%
      %=                        Find axis limits                         =%
      %===================================================================%
      zbot=min([zbot nanmin(-CTD.cruisez{i}{j}{:})]);% min(roms.cruisez{i}{j})]);
      ztop=max([ztop nanmax(-CTD.cruisez{i}{j}{:})]);% max(roms.cruisez{i}{j})]);
      smin=min([smin nanmin(CTD.cruisesalt{i}{j}) nanmin(roms.cruisesalt{i}{j})]);
      smax=max([smax nanmax(CTD.cruisesalt{i}{j}) nanmax(roms.cruisesalt{i}{j})]);
      tmin=min([tmin nanmin(CTD.cruisetemp{i}{j}) nanmin(roms.cruisetemp{i}{j})]);
      tmax=max([tmax nanmax(CTD.cruisetemp{i}{j}) nanmax(roms.cruisetemp{i}{j})]);
    
    end
    %=====================================================================%
    %=               Adjust observational temperature plot               =%
    %=====================================================================%
    axes(ax1);
    caxis([tmin tmax]);
    datetick;
    xlims=get(gca,'xlim');
    xlims=[min([CTD.cruisetime(i,1:length(CTD.cruisetemp{i})) ...
               roms.cruisetime(i,1:length(CTD.cruisetemp{i}))]) ...
           max([CTD.cruisetime(i,1:length(CTD.cruisetemp{i})) ...
               roms.cruisetime(i,1:length(CTD.cruisetemp{i}))])];
    xdist=xlims(2)-xlims(1);
    axis([xlims(1)-xdist/30 xlims(2)+xdist/30 zbot-5 ztop+1]);
    box on; 
    set(ax1,'linewidth',1,'fontsize',8,'xticklabel','');
    TobsText=text(xlims(1),zbot-zbot/14,['CTD : ',cruisenames{i}(1:10)],...
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
    %=                 Adjust observational salinity plot                =%
    %=====================================================================%
    axes(ax3);
    if smin<0; smin=2; end; 
    caxis([smin smax]);
    datetick;
    axis([xlims(1)-xdist/30 xlims(2)+xdist/30 zbot-5 ztop+1]);
    box on; 
    set(ax3,'linewidth',1,'fontsize',8,'xticklabel','');
    SobsText=text(xlims(1),zbot-zbot/14,['CTD : ',cruisenames{i}(1:10)],...
                  'HorizontalAlignment','left',...
		          'fontsize',10,'fontweight','bold');

    %=====================================================================%
    %=                    Adjust model salinity plot                     =%
    %=====================================================================%
    axes(ax4);
    caxis([smin smax]);
    datetick;
    axis([xlims(1)-xdist/30 xlims(2)+xdist/30 zbot-5 ztop+1]);
    box on; 
    set(ax4,'linewidth',1,'fontsize',8);
    SmodText=text(xlims(1),zbot-zbot/14,'US East ROMS',...
                  'HorizontalAlignment','left',...
		          'fontsize',10,'fontweight','bold');
		  
    %=====================================================================%
    %=                     Make temperature colorbar                     =%
    %=====================================================================%
    axes(cbAx1);
    cb1=colorbar('position',[0.86 0.545 0.04 0.445]);
    caxis([tmin tmax]);
    set(cb1,'linewidth',1,'fontsize',8);
    cbText1=text(2,.5,'Temperature ( ^{o}C)','fontsize',12,'rotation',270,...
                    'HorizontalAlignment','center');
                
    %=====================================================================%
    %=                     Make salinity colorbar                        =%
    %=====================================================================%
    axes(cbAx2);
    cb2=colorbar('position',[0.86 0.04 0.04 0.445]);
    caxis([smin smax]);
    set(cb2,'linewidth',1,'fontsize',8);
    cbText2=text(2,.5,'Salinity (PSU)','fontsize',12,'rotation',270,...
                     'HorizontalAlignment','center');

    %=====================================================================%
    %=                     Export figure to file                         =%
    %=====================================================================%
    set(gcf,'color','w');
    wysiwyg;
    eval(['export_fig /home/actodd/PLOTS/USeast/validate/CTD/',...
           cruisenames{i}(1:8),'.png -png -r150 -painters']);
    
    %=====================================================================%
    %=          Make figure with cruise track and model points           =%
    %=====================================================================%
    figure(2); clf;
    
    %--- Contour land mask from model and model depths ---%
    contour(roms.lon,roms.lat,roms.mask,[0 0],'k','linewidth',2); hold on
    contour(roms.lon,roms.lat,roms.h,[50 100 250 500 7500 1000:1000:4000],...
            'color',[.4 .4 .4],'linewidth',0.75);

    %--- Plot ROMS points as plus signs ---%
    plot(roms.cruiselon( i,1:length(CTD.cruisetemp{i})),...
         roms.cruiselat( i,1:length(CTD.cruisetemp{i})),'m+');
   
    %--- Plot cruise track and CTD cast sites ---%
    plot(CTD.cruiselon( i,1:length(CTD.cruisetemp{i})),...
         CTD.cruiselat( i,1:length(CTD.cruisetemp{i})),'k','linewidth',1);
    scatter(CTD.cruiselon( i,1:length(CTD.cruisetemp{i})),...
            CTD.cruiselat( i,1:length(CTD.cruisetemp{i})),25,...
            CTD.cruisetime(i,1:length(CTD.cruisetemp{i})),'filled');
    caxis([CTD.cruisetime(i,1) CTD.cruisetime(i,length(CTD.cruisetemp{i}))])  
    
    %--- Add colorbar for time scale --%
    hCbar=colorbar;
    if (diff(CTD.cruisetime(i,[1 length(CTD.cruisetemp{i})]))>2)
      cbartickstyle=6;
    else
      cbartickstyle=15;    
    end
    set(hCbar,'linewidth',1,'fontsize',8,'yticklabel',...
        datestr(get(hCbar,'ytick'),cbartickstyle));
    
    %--- Set axis and figure properties ---%
    axis xy; axis image;
    axis([min(CTD.cruiselon(i,1:length(CTD.cruisetemp{i})))-1 ...
          max(CTD.cruiselon(i,1:length(CTD.cruisetemp{i})))+1 ...
          min(CTD.cruiselat(i,1:length(CTD.cruisetemp{i})))-1 ...
          max(CTD.cruiselat(i,1:length(CTD.cruisetemp{i})))+1])
    
    set(gcf,'color','w');
    set(gca,'linewidth',1,'fontsize',10,'fontname','Helvetica');
    
    %--- Set title ---%
    title(['Cruise ',cruisenames{i}(1:8)],'fontsize',12,'fontname',...
          'Helvetica','fontweight','bold');
      
    %--- Export figure ---%
    wysiwyg;
    eval(['export_fig /home/actodd/PLOTS/USeast/validate/CTD/',...
           cruisenames{i}(1:8),'_track.png -png -r150 -painters']);

  end
end
    
