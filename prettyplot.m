function prettyplot(name,xlab,ylab)
%Makes plots look less shitty
    ax=gca;
    fig=gcf;
    title(name,'Interpreter','Latex');
    xlabel(xlab,'Interpreter','Latex');
    ylabel(ylab,'Interpreter','Latex');
    set(ax,'FontSize',12);
    set(ax,'TitleFontSizeMultiplier',2);
    set(ax,'LabelFontSizeMultiplier',1.5);
    set(ax,'TickLabelInterpreter','Latex');
    set(fig,'Color',[1 1 1]);
    colororder(ax,{'#000000','#a30f0f','#261eb3','#b59826','#7D7D7D'});
    set(fig,'DefaultLineLineWidth',2);
    lines = findobj(ax, 'Type', 'line');% Via Matlab answers user "Image Analyst"
    L=size(lines,1);
    for i=1:L
        lines(i).LineWidth=2;
    end
    
end