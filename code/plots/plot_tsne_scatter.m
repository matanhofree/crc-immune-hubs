function zaxes = plot_tsne_scatter(ydata,clustC,colorSet,sliceVal,inOpts)

    if ~exist('colorSet','var')
        colorSet = [];
    end

    defaultOpts.doPosTxt = 0;
    defaultOpts.pSize = 20;
    defaultOpts.doSort = 1;
    defaultOpts.newPlot = 1;
    defaultOpts.isCategory = 1;
    defaultOpts.yAxis_unit = ' ';
    defaultOpts.cmap = 'jet';
    defaultOpts.fontSize = 30;
    defaultOpts.symb = [];
    defaultOpts.grey =  [0.6627 0.6627 0.6627];
    % defaultOpts.doAlpha = 0.7;
    defaultOpts.doAlpha = [];
    
    defaultOpts.includeNonSlice = 1;
    defaultOpts.nonSliceTxt = 'Other';
    defaultOpts.doLeg = 'on';
    
    defaultOpts.cmapLimits = [-5 5];
    defaultOpts.cmapQuant = [ 0.01 0.99 ];
    defaultOpts.cmapLimitsRemoveNeg = 1;
    defaultOpts.plotSize = [ 10 10 1200 1100];
    defaultOpts.useGramm= 0;
    defaultOpts.dotText = [];
    defaultOpts.dotTextSize = 9;
            
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    
    if opts.cmapLimitsRemoveNeg && ~iscell(clustC) && all(clustC >= 0)
        opts.cmapLimits = max(opts.cmapLimits,0);
    end
        
    % disp(opts);
    if opts.newPlot 
        zfig = figure('Position',opts.plotSize );
        set(gcf,'color','w'); 
    end
    zaxes = gca();
                   
    clustC = clustC(:);
    
    % If sorting just sort
    if opts.doSort == 1
        [~,zOidx] = sort(clustC);
        
        clustC = clustC(zOidx);
        ydata = ydata(zOidx,:);
    elseif length(opts.doSort)>1
        zOidx = opts.doSort;
        
        clustC = clustC(zOidx);
        ydata = ydata(zOidx,:);
    else
        zOidx = 1:length(clustC);
    end  
    
    
    
    sliceSelect = [];
    if exist('sliceVal','var') && ~isempty(sliceVal)
        if length(sliceVal) == length(clustC)                        
            sliceSelect = sliceVal(zOidx);
        elseif ischar(sliceVal)
            sliceSelect = strgrep(clustC,sliceVal);     
        else
            error('Size missmatch');
        end            
        
        % ydataOther = ydata(~sliceSelect,:);
        % clustCother = clustC(~sliceSelect);
    end
      

    if ~isempty(sliceSelect)
        
        colorOrder = uniquec(clustC);    
        % zOtherColor = unique(clustCother);
        zOtherColor = colorOrder;
        zOtherColor(ismember(zOtherColor,unique(clustC(sliceSelect)))) = [];
    end 
    
    
    if opts.isCategory && (iscell(clustC) || length(unique(clustC)) < 100)
        
        if isempty(colorSet)
            colorSet = jet(luniq(clustC));            
        elseif isa(colorSet,'containers.Map')
            colorSet = nanvalues(colorSet,colorOrder);
            
            colorSet(isemptycell(colorSet)) =  {opts.grey};
            
            colorSet = cell2mat(colorSet');
        else            
            if ~isempty(sliceSelect)
                % colorOrder = uniquec(clustC);    

                sel = ismember(colorOrder,zOtherColor);
                colorSet(sel,:) = repmat(opts.grey,sum(sel),1);
            end
        end
        
        if isempty(sliceSelect)            
            if opts.useGramm
                outG = gramm_gscatter(ydata(:,1),ydata(:,2),clustC,colorSet,opts.symb,opts.pSize,opts.doLeg);
            else
                gscatter(ydata(:,1),ydata(:,2),clustC,colorSet,opts.symb,opts.pSize,opts.doLeg);
                % zaxes.Parent.Children(1).Interpreter = 'none';
                zf = legend('boxoff');
                zf.Interpreter = 'None';
                
            end
            
            
        else

            
            ydataDrop = ydata;
            ydataDrop(sliceSelect,:) = nan;                        
            if ~isempty(opts.doAlpha) && opts.doAlpha > 0
                colorSetAlpha = brighten(colorSet,opts.doAlpha);
                % colorSetAlpha(:,4) = opts.doAlpha;
                zfig = gscatter(ydataDrop(:,1),ydataDrop(:,2),clustC,colorSetAlpha,opts.symb,opts.pSize,'off');
            else
                zfig = gscatter(ydataDrop(:,1),ydataDrop(:,2),clustC,opts.grey,opts.symb,opts.pSize,'off');
            end
            
            hold on;
            ydataDrop = ydata;
            ydataDrop(~sliceSelect,:) = nan;
            gscatter(ydataDrop(:,1),ydataDrop(:,2),clustC,colorSet,opts.symb,opts.pSize,opts.doLeg);
            % zaxes.Parent.Children(1).Interpreter = 'none';
            legend('boxoff');

%             zz.Children(2).Children = flipud(zz.Children(2).Children);
%             zz.Children(2).Children = [ flipud(zz.Children(2).Children(~sel)); zz.Children(2).Children(sel);];
            
        end 
        
        
    else
        if ~isempty(opts.cmapLimits) || ~isempty(opts.cmapQuant)
            
            if ~isempty(opts.cmapQuant)
                if all(opts.cmapLimits<1) && all(opts.cmapLimits>=0) && length(opts.cmapQuant) == 1
                    opts.cmapLimits = quantile(clustC,opts.cmapLimits);               
                elseif length(opts.cmapQuant) == 2
                    opts.cmapLimits = quantile(clustC,opts.cmapQuant);               
                end
            end
            
            clustCmaxed = max(min(clustC,opts.cmapLimits(2)),opts.cmapLimits(1));

        else
            clustCmaxed = clustC;
        end
        
        if ~isempty(sliceSelect)                
            zf = scatter(ydata(~sliceSelect,1),ydata(~sliceSelect,2),opts.pSize,opts.grey,'filled');
            hold on;
            zf = scatter(ydata(sliceSelect,1),ydata(sliceSelect,2),opts.pSize,clustCmaxed(sliceSelect),'filled');
        else            
            zf = scatter(ydata(:,1),ydata(:,2),opts.pSize,clustCmaxed,'filled');
        end
        
        axis tight
        clim = ylim();
        ylim([ clim(1) - diff(clim)*0.025 clim(2) + diff(clim)*0.025])
        clim = xlim();
        xlim([ clim(1) - diff(clim)*0.025 clim(2) + diff(clim)*0.025])
        
        colormap(opts.cmap); 
        if ~isempty(opts.cmapLimits) && ~all(opts.cmapLimits ==0)
            % c.Limits = opts.cmapLimits;
            c = colorbar('Limits',opts.cmapLimits); 
            % caxis(opts.cmapLimits);
        else
            c = colorbar(); 
        end
        
        ylabel(c,opts.yAxis_unit,'interpreter','none');
        
    end
    
    xlabel('tSNE 1'); ylabel('tSNE 2');
    box off 
    
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);

    if (opts.doPosTxt)
        clustList = unique(clustC);
        
        if iscell(clustList)
            
            meanCt_x = cellfun(@(x)median(ydata(strcmp(clustC,x),1)),clustList);
            meanCt_y = cellfun(@(x)median(ydata(strcmp(clustC,x),2)),clustList);
            text(meanCt_x,meanCt_y,clustList,'fontsize',opts.fontSize,'HorizontalAlignment','center','interpreter','none');
        else
            meanCt_x = arrayfun(@(x)median(ydata(clustC == x,1)),clustList);
            meanCt_y = arrayfun(@(x)median(ydata(clustC == x,2)),clustList);
            text(meanCt_x,meanCt_y,num2cellstr(clustList),'fontsize',opts.fontSize,'HorizontalAlignment','center');
        end
    end
    
    if ~isempty(opts.dotText)
        if length(opts.dotText) ~= length(ydata)
            fprintf('Dot annotation size missmatch');
        else
            cText = opts.dotText(zOidx);
            text(ydata(:,1),ydata(:,2),cText);
            
        end
    end

    plot_dump(mfilename());    
end