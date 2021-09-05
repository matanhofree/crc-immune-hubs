function gg = plot_graph_coAct(G,cNames,cNodesSize,cEdges,cMap,cMeanE,cMeanStd,cFDR)

    zsG = @(x,me,st)((x-me)/st);

    zTrimNames = @(x)regexprep(x,'.*_zz_','');
    
    
    % cLWidths = max(3 + 5*zsG(cEdges,cMeanE,cMeanStd),1)
    cLWidths = max(30*abs(cEdges),3);

    gg = plot(G,'Interpreter','None','layout','force','usegravity','on','Iteration',10000,'NodeLabel',zTrimNames(cNames),'NodeColor',cMap,'MarkerSize',cNodesSize,'LineWidth',cLWidths,'WeightEffect','direct');

    
    zCmap = [    0.9569    0.5843    0.7490;
                0.1137    0.4627    0.7333];

    zCol = zeros(length(cEdges),3);
    zCol(cEdges < 0,:) = repmat(zCmap(end,:),sum(cEdges < 0),1);
    zCol(cEdges > 0,:) = repmat(zCmap(1,:),sum(cEdges > 0),1);

    gg.EdgeColor = zCol;
    gg.NodeFontSize = 16;


    zSel = cFDR < 0.1;
    
    zLineSt = repmat({':'},1,length(zSel));
    zLineSt(zSel) = { '-' };
    
    gg.LineStyle = zLineSt;
    
end

