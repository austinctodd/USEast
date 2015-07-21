a=rand(600,400);
b=rand(600,400);
c=b-a;

figure(1); clf

%--- Set axes and plot frame 1 ---%
%position in normalized size (0 to 1)
%bracket is [X_start,Y_start,Width,Height]
h1=axes('position',[0.12 0.5 0.425 0.425]);
imagesc(a'); axis image;
set(gca,'fontsize',12,...              % set font size of tick labels
        'xtick',[100:100:500],...      % set x tick marks
        'xticklabel','',...            % remove x tick mark labels
        'ytick',[50:50:350]);          % set y tick marks
ylabel('Pixels','FontSize',14);
    
%--- Set axes and plot frame 2 ---%
%position in normalized size (0 to 1)
%bracket is [X_start,Y_start,Width,Height]
h2=axes('position',[0.56 0.5 0.425 0.425]);
imagesc(b'); axis image;
set(gca,'fontsize',12,...              % set font size of tick labels
        'xtick',[100:100:500],...      % set x tick marks
        'xticklabel','',...            % remove x tick mark labels
        'ytick',[50:50:400],...        % set y tick marks
        'yticklabel','');              % remove y tick mark labels

%--- Set axes and plot difference ---%
%position in normalized size (0 to 1)
%bracket is [X_start,Y_start,Width,Height]
h3=axes('position',[0.12 0.1 0.425 0.425]);
imagesc(c'); axis image;
set(gca,'fontsize',12,...              % set font size of tick labels
        'xtick',[100:100:500],...      % set x tick marks
        'ytick',[50:50:350]);           % set y tick marks
xlabel('Pixels','FontSize',14);
ylabel('Pixels','FontSize',14);

%--- Set axes and plot contour contraints ---%
%position in normalized size (0 to 1)
%bracket is [X_start,Y_start,Width,Height]
h4=axes('position',[0.61 0.127 0.4 0.37]);
imagesc(415:455,215:247,c(415:455,215:247)'); axis image;
set(gca,'fontsize',12,...              % set font size of tick labels
        'xtick',[420:10:450],...       % set x tick marks
        'ytick',[220:5:245]);           % set y tick marks
xlabel('Pixels','FontSize',14);
ylabel('Pixels','FontSize',14);


