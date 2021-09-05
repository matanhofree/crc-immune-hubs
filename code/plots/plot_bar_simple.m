function [zfig,outG,g] = plot_bar_simple(cXi,cYi,cFacti,inOpts)
    
    defaultOpts.newPlot = 1;
    defaultOpts.maxRow = 20;
    defaultOpts.widthV = 0.6;
%     defaultOpts.widthBox = 0.3;
    defaultOpts.dodge = 0.4;
%    defaultOpts.outliers = 0;
%    defaultOpts.whiskers = 0;
%    defaultOpts.notch = 0;
    defaultOpts.doJitter = 0;
    defaultOpts.ylim = [];
    defaultOpts.plotSize = [1 1 1600 900];
    defaultOpts.reorder = 0;
    defaultOpts.reorderX = 0;
    defaultOpts.cmap = [];
    defaultOpts.reorderXbyY = 0;
    
    defaultOpts.title = [];
    defaultOpts.groupMode = 1;
    defaultOpts.dropLegend = 0;
    defaultOpts.yLabel = [];
    defaultOpts.xLabel = [];
    defaultOpts.font = 'Ariel';
    defaultOpts.fontSize = 16;
    
    defaultOpts.XTickLabelRotation = [];
    
    defaultOpts.errorBar = [];
    
    defaultOpts.flipCoor = 0;
    defaultOpts.drawFig = 1;
    
    zfig = [];
    outG = [];
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts)
     
        
    if ~exist('cYi','var') || isempty(cYi)
        cYi = [];               
    end
    if ~exist('cFacti','var') || isempty(cFacti)
        cFacti = [];               
    end          
    
    zfig = [];
    if opts.newPlot == 1
        zfig = figure('position',opts.plotSize);        
    end   
    
    if istable(cXi)   
        
        cY = table2array(cXi);        
        cFact = cXi.Properties.VariableNames;
        cX = cXi.Properties.RowNames;
       
        [N,nC] = size(cY);
        %%
        
        if ~isempty(cFacti) && length(cFacti) == nC
            cFact = cFacti;
        end
            
            
        [cFact,cX] = meshgrid(cFact,cX);
        
        if ~isempty(cFacti) && isnumeric(cFacti) && all(cFacti(:) == 0)
            cFact = []
        end
        
            
    else        
%         [N,nC] = size(cXi);
%         if nargin<3 || isempty(cFacti)
%             cFact = [];                        
%         else
%             cFact = matchVector(cXi,cFacti,N,nC);
%         end


        % cY = matchVector(cX,cY,N,nC);            
        
        if nargin<3 || isempty(cFacti)
            cFact = [];                        
        else
            cFact = matchVector(cYi,cFacti);
        end
        
        cX = matchVector(cYi,cXi);             
        cY = cYi;
    end
    
    
        
    cError = [];
    if ~isempty(opts.errorBar)
       cError = [ cY-opts.errorBar cY+opts.errorBar ];
    end
    
        
    cX = cX(:);
    cY = cY(:);

    if ~isempty(cFact)
        cFact = cFact(:);
    end
    
    [cXlist,~,~,~,cXpos] = fastUnique(cX);
    if opts.reorderX
        [cXlist,xReord] = sort(cXlist);
        
        cXpos = cXpos(xReord);        
    end
    
    nC = length(cXlist);
    nR = ceil(nC/opts.maxRow); 
    
    zp = 1;    
    for i = 1:nR
        
        cP = zp:min(zp+opts.maxRow-1,nC);
        rowLen = length(cP);
        
        cSub = [ cXpos{cP} ];   
        
        if isempty(cFact)
            g(i,1) = gramm('x',cX(cSub),'y',cY(cSub));
        else
            if ~isempty(cError)
                g(i,1) = gramm('x',cX(cSub),'y',cY(cSub),'color',cFact(cSub),'ymin',cError(cSub,1),'ymax',cError(cSub,2));
            else
                g(i,1) = gramm('x',cX(cSub),'y',cY(cSub),'color',cFact(cSub));
            end
        end
        
        if ~isempty(opts.widthV)
            g(i,1).stat_summary('type','sem','geom',{'bar','black_errorbar'},'width',opts.widthV,'dodge',opts.dodge);                        
        end
        
        if ~isempty(cError)
            g(i,1).geom_interval('geom','black_errorbar','width',opts.widthV,'dodge',opts.dodge);
        end
%         if ~isempty(opts.widthBox)
%             g(i,1).stat_boxplot('width',opts.widthBox,'dodge',opts.dodge,'outliers',opts.outliers,'whisk',opts.whiskers,'notch',opts.notch);
%         end

        if opts.doJitter > 0
            g(i,1).geom_jitter('dodge',opts.dodge,'width',opts.doJitter);
        end
        if ~isempty(opts.ylim)
            g(i,1).axe_property('YLim',opts.ylim);
        end
        if opts.reorder == 0
            g(i,1).set_order_options('x',0);
        elseif opts.reorderXbyY 
            error('Not implemented');
%             if ~isempty(cFact)
%                 [zFactList,~,~,~,zFactOrder] = fastUnique(cFact(cSub));
%                 
%             else
%                 [zMlist,~,~,~,zMpos] = fastUnique();            
%                 cMedianY =   
%             emd
            
        end        
                  
        g(i,1).set_text_options('base_size',12);
        if ~isempty(opts.title) && i == 1;
            g(i,1).set_title(opts.title);            
        end
        
        if ~isempty(opts.yLabel) || strcmp(opts.yLabel,'')
            g(i,1).set_names('y',opts.yLabel);            
        end
        
        if ~isempty(opts.xLabel) || strcmp(opts.xLabel,'')
            g(i,1).set_names('x',opts.xLabel);            
        end
        
        
        zp = zp + opts.maxRow;
    end
    if ~exist('g','var') || isempty(g)
        warning('Did not create a plot!');
        g = [];
        return;
    end
        
    if opts.dropLegend == 1
        g.set_layout_options('legend',false);
    end
    if ~isempty(opts.cmap)
        g.set_color_options('map',opts.cmap,'n_color',size(opts.cmap,1));
    end    
    g.set_text_options('font',opts.font,'base_size',opts.fontSize);
    
    if opts.flipCoor
        g.coord_flip();
    end
    if opts.drawFig
        outG = g.draw();
    end
    
    if ~isempty(opts.XTickLabelRotation)
        for zzc = 1:length(outG)
            outG(zzc).facet_axes_handles.XTickLabelRotation = opts.XTickLabelRotation;
        end
        outG = outG.redraw();
    end
    
    
    
end


function cY = matchVector(refX,eV)    
    
    if ~all(size(refX) == size(eV))
        
        [nr,nc] = size(refX);
        if (nr == nc) 
            fprintf('Dimensions are equal dims must be expanded manually');
            error('Unable to expand matching dims');
        end
        
        if min(size(eV)) ~= 1
            error('Unable to expand non vector dims');
        end
        
        if  nr == length(eV)
            cY = repmat(eV(:),1,nc);
        elseif nc == length(eV)
            cY = repmat(eV(:)',nr,1);
        else
            error('Size missmatch on %s',inputname(1));            
        end            
    else
        cY = eV;
    end

end
