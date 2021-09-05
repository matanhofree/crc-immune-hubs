function [figOut,outMat,zOrdTX,zOrdTY] = plot_heatmap_annot(crossTable,yName,xName,cmap,textTable,inOpts)

    defaultOpts.axisObj =[];
    defaultOpts.showTxt = 0;
    defaultOpts.XTickLabelRotation = 60;

    defaultOpts.doSortY = 1;
    defaultOpts.doSortX = 1;
    defaultOpts.doLeafOptimalOrderY = 2;
    defaultOpts.doLeafOptimalOrderY_data = [];
    
    defaultOpts.doLeafOptimalOrderX = 2;
    defaultOpts.doLeafOptimalOrderX_data = []; 
        
    defaultOpts.externXOrder = 0;
    defaultOpts.externYOrder = 0;
    
    defaultOpts.colormap = redblue(128); % flipud(cptcmap('temperature')); % flipud(hsv(128)); % redblue(128); % redbluecmap();% flipud(cptcmap('ylwhbl_sym'));
    defaultOpts.colorLim = [ 0.01 0.99 ];
    defaultOpts.colorHard = [];
    
    defaultOpts.fontSize = 12;
    defaultOpts.axisFontSize = 12;
    
    % defaultOpts.leafDistTrans = 'linear';
    defaultOpts.optLeafDistTrans = 'linear';
    defaultOpts.optLeafCrit = 'adjacent';
    defaultOpts.doOptLeaf = 1;
    
    defaultOpts.addDendrogramX = 1; 
    defaultOpts.addDendrogramY = 1;
    
    defaultOpts.linkageX = 'average';
    defaultOpts.linkageY = 'average';
    
    defaultOpts.pdistX = 'corr';
    defaultOpts.pdistY = 'corr';
    
    defaultOpts.pseudoCount = 1;

    defaultOpts.maxLabelX = 100;
    defaultOpts.maxLabelY = 100;
    
    defaultOpts.outClust = 0;

    defaultOpts.annotX = [];
    defaultOpts.annotY = [];
    defaultOpts.textPrec = 2;
    
    defaultOpts.colorSym = 1;
    
    defaultOpts.position = [ 1 1 1800 1800 ];
    
    defaultOpts.subAxisOpts = { 'SpacingVertical',0.01,'SpacingHorizontal',0.01 };

    
    defaultOpts.useUpperLeftAnnot = 1;
    
    defaultOpts.blockSortY = 0;
    defaultOpts.doDotPlot = 0;
    defaultOpts.dotPlotFactor = 1000;
    defaultOpts.dotPlotColor = [ 0 0 0];
    defaultOpts.dotPlotLineWidth = 1;
    
    defaultOpts.leafOptimalOrderX_const = [];
    defaultOpts.gridLines = 'none';
%     pLeft = 0.01;
%     pBottom = 0.05;
    
                
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    % disp(opts);
    zfig = [];
    
    if istable(crossTable)
        if ~exist('yName','var') || isempty(yName)
            yName = crossTable.Properties.RowNames;
        end
        
        if ~exist('xName','var') || isempty(xName)
            xName = crossTable.Properties.VariableNames;
        end
        
        crossTableO = crossTable;
        crossTable = table2array(crossTable);
    end
    
    
    if ~exist('textTable','var') || isempty(textTable)
        textTable = crossTable;        
    end
    
    outMat = [];
    zOrdTX = []
    zOrdTY = [];
          
    if isempty(opts.axisObj)                     
        figOut = figure('color','w','Position', opts.position); 
        
        if opts.addDendrogramX > 0 || opts.addDendrogramY > 0  
            
%             if isempty(opts.annotX) && isempty(opts.annotY)  
%                 
%                 % dendyP = subaxis(4,4,1,2,1,3);
%                 % dendxP = subaxis(4,4,2,1,3,1);      
%                 heatMapP = subaxis(4,4,2,2,3,3);                                                  
%                 
%             else        
                if opts.useUpperLeftAnnot == 1
                    heatMapLegendSpace = subaxis(12,12,1,1,3,3,opts.subAxisOpts{:});
                   
%                 elseif opts.useUpperLeftAnnot == 2
%                     
%                 elseif opts.useUpperLeftAnnot == 3
%                     
                    
                end
                
                dendyP = subaxis(12,12,1,4,2,9,opts.subAxisOpts{:});
                dendxP = subaxis(12,12,4,1,9,2,opts.subAxisOpts{:}); 
                
                if ~isempty(opts.annotX)                    
                    posAnnotX = subaxis(12,12,4,3,9,1,opts.subAxisOpts{:}); 
                end
                
                if ~isempty(opts.annotY)
                    posAnnotY = subaxis(12,12,3,4,1,9,opts.subAxisOpts{:}); 
                end
                
                heatMapP = subaxis(12,12,4,4,9,9,opts.subAxisOpts{:});  
                
                
            % end
                        
        else            
            heatMapP = gca;
        end
    elseif isstruct(opts.axisObj)
        
        if isfield(opts.axisObj,'heatmap')
            heatMapP = opts.axisObj.heatmap;
        end
        if isfield(opts.axisObj,'treeY')
            dendyP = opts.axisObj.treeY;
        end
        if isfield(opts.axisObj,'treeX')
            dendxP = opts.axisObj.treeX;
        end
        
        
    else
        heatMapP = opts.axisObj;
    end         
    
    [nY,nX] = size(crossTable);
        
    if ~exist('yName','var') || isempty(yName)        
        yName = 1:nY;
        opts.doSortY = 0;
    end

    if ~exist('xName','var') || isempty(xName)        
        xName = 1:nX;
        opts.doSortX = 0;
    end

    outMat.xName = xName;
    outMat.yName = yName;
    
    if ~exist('cmap','var') || isempty(cmap)        
        cmap = opts.colormap;
        colormap(cmap);
    end
    
    xOrder = [];
    yOrder = [];
    leafOrderY = [];
    leafOrderX = [];
    if opts.doSortY 
        if isnumeric(yName)
            [~,yOrder] = sort(str2double(yName));
        else
            [~,yOrder] = sort(yName);
        end                
        yName = yName(yOrder);
        crossTable = crossTable(yOrder,:);
        textTable = textTable(yOrder,:);
    elseif opts.externYOrder 
        fprintf('Using external Y order\n');
        yOrder = opts.externYOrder;
        yName = yName(yOrder);
        crossTable = crossTable(yOrder,:);
        textTable = textTable(yOrder,:);
    end
    
    if opts.doSortX
        if isnumeric(xName)
            [xName,xOrder] = sort(str2double(xName));
        else            
            [xName,xOrder] = sort(xName);
        end
        
        crossTable = crossTable(:,xOrder);
        textTable = textTable(:,xOrder);
    elseif opts.externXOrder 
        fprintf('Using external X order\n');
        xOrder = opts.externXOrder;
        xName = xName(xOrder);
        crossTable = crossTable(:,xOrder);
        textTable = textTable(:,xOrder);
    
    end

    %%
    % Alternative sort method 
    if opts.doLeafOptimalOrderY
        if isempty(opts.doLeafOptimalOrderY_data)
            if opts.doLeafOptimalOrderY == 1 
                fprintf('Calculating similarity in freq. Y (doLeafOptimalOrderY == 1)\n');
                zDY = pdist(crossTable./sum(crossTable,2),opts.pdistY);
                zOrdTY = linkage(zDY,opts.linkageY);
            elseif opts.doLeafOptimalOrderY == 2 
                zDY = pdist(crossTable,opts.pdistY);
                zOrdTY = linkage(zDY,opts.linkageY);            
            end
        else            
            zDY = opts.doLeafOptimalOrderY_data.pdist;
            if size(zDY,1) ~= size(zDY,2)
                zDY = squareform(zDY);
            end

            if ~isempty(yOrder)
                zDY = zDY(yOrder,yOrder);
            end
            zOrdTY = linkage(squareform(zDY),opts.linkageY);   

        end
        
        if opts.doOptLeaf
            leafOrderY = optimalleaforder(zOrdTY,zDY,'Transformation',opts.optLeafDistTrans,'criteria',opts.optLeafCrit);
        else
            fprintf('Denrogram Order')
            [~,leafOrderY] = dendrogram_order(zOrdTY);
        end
        crossTable = crossTable(leafOrderY,:);
        textTable = textTable(leafOrderY,:);        
        yNameO = yName(leafOrderY);               
    else
        yNameO = yName;
    end
        
    if opts.doLeafOptimalOrderX
        if isempty(opts.doLeafOptimalOrderX_data)
%             if opts.doLeafOptimalOrderX == 1 
%                 fprintf('Calculating similarity in freq. X (doLeafOptimalOrderX == 1)\n');
%                 
%                 zDX = pdist((crossTable./sum(crossTable,1))',opts.pdistX);
%                 zOrdTX = linkage(zDX,opts.linkageX);
%             elseif opts.doLeafOptimalOrderX == 2 
%                 zDX = pdist(crossTable',opts.pdistX);
%                 zOrdTX = linkage(zDX,opts.linkageX);            
%             elseif opts.doLeafOptimalOrderX == 3 
%                 zDX = pdist(zscore(crossTable)',opts.pdistX);
%                 zOrdTX = linkage(zDX,opts.linkageX);            
%             elseif opts.doLeafOptimalOrderX == 4
%                 zDX = pdist(zscore(crossTable'),opts.pdistX);
%                 zOrdTX = linkage(zDX,opts.linkageX);                        
%             end
            if opts.doLeafOptimalOrderX == 1 
                fprintf('Calculating similarity in freq. X (doLeafOptimalOrderX == 1)\n');
                cD = (crossTable./sum(crossTable,1))';
            elseif opts.doLeafOptimalOrderX == 2 
                cD = crossTable';                
            elseif opts.doLeafOptimalOrderX == 3 
                cD = zscore(crossTable)';                
            elseif opts.doLeafOptimalOrderX == 4
                cD = zscore(crossTable');                            
            end

            if isempty(opts.leafOptimalOrderX_const)
                zDX = pdist(cD,opts.pdistX);
                zOrdTX = linkage(zDX,opts.linkageX);            
            else
                opts.doOptLeaf = 0;
                leafOrderX = constrainedLeafOrder(cD,opts.leafOptimalOrderX_const(xOrder),2,opts.pdistX,opts.linkageX,opts.doOptLeaf);                
            end
        else                        
            zDX = opts.doLeafOptimalOrderX_data.pdist;
            if size(zDX,1) ~= size(zDX,2)
                zDX = squareform(zDX);
            end
            
            if ~isempty(xOrder)
                zDX = zDX(xOrder,xOrder);
            end

            zOrdTX = linkage(squareform(zDX),opts.linkageX);     
        end

        if opts.doOptLeaf
            fprintf('Optimal leaf orderX\n')
            leafOrderX = optimalleaforder(zOrdTX,zDX,'Transformation',opts.optLeafDistTrans,'criteria',opts.optLeafCrit);
        elseif isempty(opts.leafOptimalOrderX_const)
            fprintf('Denrogram OrderX\n')
            [~,leafOrderX] = dendrogram_order(zOrdTX);

        end
        
        crossTable = crossTable(:,leafOrderX);
        textTable = textTable(:,leafOrderX);
        xNameO = xName(leafOrderX);
                
    else
        xNameO = xName;
    end
    
    
    if opts.blockSortY 
        if opts.doLeafOptimalOrderY
            warning('Cant block sort and draw dendrogram');
        else
        
            fprintf('Block sorting Y');
        
            [zz,mxidx] = sort(crossTable,2,'Descend');                      
            
            % [~,leafOrderY] = sortrows(mxidx);            
            % Todo add sub optimal leaf order
            opts.sortClustNames = 1;
            leafOrderY = constrainedLeafOrder(crossTable,mxidx(:,1),1,opts.pdistY,opts.linkageX,1,opts)     
            
            crossTable = crossTable(leafOrderY,:);
            textTable = textTable(leafOrderY,:);        
            yNameO = yNameO(leafOrderY)
            
        end
    end
    
    if ~isempty(opts.colorLim)
        
        clim = [ 0 0 ];
        
        if opts.colorLim(2) >= 1
            clim(2) = opts.colorLim(2);
        else
            clim(2) = quantile((crossTable(:)),opts.colorLim(2));
        end
        
        if opts.colorLim(1) > 1 || opts.colorLim(1) <= 0
            clim(1) = opts.colorLim(1);
        else
            clim(1) = quantile((crossTable(:)),opts.colorLim(1));
        end
        
    else
        clim = max(abs(crossTable(:)));
        clim = [ -clim clim ];
    end
       
    if opts.colorSym 
        clim = max(abs(clim(:)));
        clim = [ -clim clim ];
    end
    
    if ~isempty(opts.colorHard)
        clim = opts.colorHard;
    end

   

    [nY,nX] = size(crossTable);
    
    xNamePrint = xNameO;
    yNamePrint = yNameO;
    
    if ~isempty(opts.maxLabelX) && opts.maxLabelX < nX        
        xNamePrint = [];
    end
        
    if ~isempty(opts.maxLabelY) && opts.maxLabelY < nY        
        yNamePrint = [];                
    end

    if opts.showTxt == 1       
        [hImage, hText, hXText] = heatmapTXT(crossTable,xNamePrint, yNamePrint, round(textTable,opts.textPrec), 'TickAngle', opts.XTickLabelRotation,...
            'ShowAllTicks', true, 'TickFontSize', opts.axisFontSize,'FontSize',opts.fontSize,'Colormap',cmap,'Colorbar',{'location','eastoutside'},'mincolorvalue',clim(1),'maxcolorvalue',clim(2),'Parent',heatMapP,'GridLines',opts.gridLines);
    elseif opts.showTxt ==2 
        [hImage, hText, hXText] = heatmapTXT(crossTable,xNamePrint, yNamePrint, textTable, 'TickAngle', opts.XTickLabelRotation,...
            'ShowAllTicks', true, 'TickFontSize', opts.axisFontSize,'FontSize',opts.fontSize,'Colormap',cmap,'Colorbar',{'location','eastoutside'},'mincolorvalue',clim(1),'maxcolorvalue',clim(2),'Parent',heatMapP,'GridLines',opts.gridLines);            
    else
        %%
        [hImage, hText, hXText] = heatmapTXT(crossTable,xNamePrint, yNamePrint, [], 'TickAngle', opts.XTickLabelRotation,...
            'ShowAllTicks', true, 'TickFontSize', opts.axisFontSize,'FontSize',opts.fontSize,'Colormap',cmap,'Colorbar',{'location','eastoutside'},'mincolorvalue',clim(1),'maxcolorvalue',clim(2),'Parent',heatMapP,'GridLines',opts.gridLines);       
    end
    %%
    cax = hImage.Parent;
    set(cax,'TickLength',[0 0])

    if opts.doDotPlot
        
        cSize = textTable';

        [zi,zj,zv] = find(cSize);
        
        zv = zv*opts.dotPlotFactor;
        %%

        cax = hImage.Parent;
        hold on 
        %
        cax2 = scatter(zi,zj,zv,'Parent',cax);
        cax2.CData = opts.dotPlotColor;
        cax2.LineWidth = opts.dotPlotLineWidth;
        
        % set(cax,'TickLength',[0 0])
    end
   
    if opts.addDendrogramY  > 0 && opts.doLeafOptimalOrderY
                
        axes(dendyP);
        zfdY = dendrogram(zOrdTY,0,'Orientation','left','Reorder',flipud(leafOrderY(:)));     
        axis off
        % axis tight
        
        dendyP.YLim = [ min(dendyP.YTick) - 0.5 max(dendyP.YTick) + 0.5 ];
    end    

    %%
    if opts.addDendrogramX > 0 && opts.doLeafOptimalOrderX && exist('zOrdTX','var') && ~isempty(zOrdTX)
        
        dendxP.Position(3) = hImage.Parent.Position(3);
        axes(dendxP);
        zfdX = dendrogram(zOrdTX,0,'Orientation','Top','Reorder',leafOrderX(:));
        axis off
        % axis tight
        dendxP.XLim = [ min(dendxP.XTick) - 0.5 max(dendxP.XTick) + 0.5 ];
    end    
    
    
    %%
    if ~isempty(opts.annotX)
        %%
        outAnnotX.sampleID = xNameO(:);
        outAnnotX = addSampleSlice(outAnnotX,opts.annotX,[],1,1);
        %%
        posAnnotX.Position(3) = hImage.Parent.Position(3);
        axes(posAnnotX);
        
        %%
        legendAxis = [];
        if opts.useUpperLeftAnnot == 1
            
            lN = length(fieldnames(outAnnotX)) - 1;
            
            if lN ==  1
                legendAxis = heatMapLegendSpace;
            elseif lN == 2
              
                cPos = heatMapLegendSpace.Position;
                xD = heatMapLegendSpace.Position(3)-heatMapLegendSpace.Position(1);
                
                cPosA = cPos;
                cPosB = cPos;
                
                cPosA(3) = xD/2 - 0.001;
                cPosB(3) = xD/2 - 0.001;
                cPosB(1) = cPos(1) + xD/2 + 0.001;
                legendAxis{1} = subplot('Position',cPosA);
                legendAxis{2} = subplot('Position',cPosB);
                
                
            elseif lN == 3
                error('Use external legend -- could still be done perhaps.');
            else
                error('Use external legend');
            end                            
        end
            
        [ax,legendFig] = plot_add_annot(posAnnotX,rmfield(outAnnotX,'sampleID'),1,legendAxis,opts)
    end
    
    
    %% Fixup outputs
    % outMat.cTable = array2table(crossTable,'VariableNames',matlab.lang.makeUniqueStrings(zRefNameO'),'rownames',matlab.lang.makeUniqueStrings(xNameO));
    % outMat.cPR = array2table(cPR,'variablenames',matlab.lang.makeUniqueStrings(zRefNameO), 'rownames',matlab.lang.makeUniqueStrings(xNameO));
    if exist('zOrdTX','var') && opts.outClust 
        outMat.distX = zDX;
        outMat.ordTX = zOrdTX;
    end
    if exist('zOrdTY','var') && opts.outClust 
        outMat.distY = zDY;
        outMat.ordTY = zOrdTY;        
    end
    outMat.outTable = crossTable;
    %if ~isequal(crossTable,crossTable)
%        outMat.cPR = crossTable;
%    end

    if ~isempty(textTable)
        outMat.textTable = textTable;
    end
%     outMat.xName = xName;
%     outMat.yName = yName;
    outMat.xNameOrd = xNameO;
    outMat.yNameOrd = yNameO;
    if ~isempty(xOrder) && ~isempty(leafOrderX)
        outMat.xOrder = xOrder(leafOrderX);
    else
        if exist('leafOrderX','var')
            outMat.xOrder = (leafOrderX);
        else
            outMat.xOrder = 1:length(xNameO);
        end
    end
    if ~isempty(yOrder)  && ~isempty(leafOrderY)
        outMat.yOrder = yOrder(leafOrderY);
    else
        if exist('leafOrderY','var')
            outMat.yOrder = (leafOrderY);
           
        else
            outMat.yOrder = 1:length(yNameO);
        end
    end

    %%
    % plot_dump([mfilename() '_main'],zfig);
    
end

