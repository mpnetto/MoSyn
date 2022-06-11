function plotWeiREA(app, AdjMat, COORD, Elect)

    weights = AdjMat(:,3);
    Edges = AdjMat(:,[1,2]);
    NormWei = weights./max(weights);

    x = COORD(:,1); y = COORD(:,2);

    cm = colormap('winter');
    if ~ishold 
        hold on
    end
    for yy=1:length(NormWei)
        colorID(yy,1) = max(1, sum(NormWei(yy,1) > [0:1/length(cm(:,1)):1]));
    end
    
    R = (max(COORD(:,1))-min(COORD(:,1)))/2+0.05; ang=0:0.01:2*pi;
    xp=R*cos(ang); yp=R*sin(ang);
    xn = [0.4 0.5 0.6]; yn = [0.5+R-0.03 0.5+R+0.05 0.5+R-0.03];
    plot( 0.5+xp,0.48+yp,'k','LineWidth',2);
    hold on
    plot(xn,yn,'k','LineWidth',2);

    myColor = cm(colorID, :);
    set(gca,'ColorOrder',myColor)
    Hline = plot(x(Edges).',y(Edges).','-');
    set(Hline, {'LineWidth'}, num2cell((NormWei)*20))
    plot(app, COORD(:,1),COORD(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor',[0,0,0]);
    
    text(COORD(:,1) -0.03, COORD(:,2) - 0.04, Elect, ...
        'FontSize', 10');

end