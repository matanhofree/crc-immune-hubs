function [zfig,outG,g] = plot_violin_simple(cX,cY,cFact,inOpts)
    
    defaultOpts.newPlot = 1;
    defaultOpts.maxRow = 20;
    defaultOpts.nSplit = 0;
    defaultOpts.widthV = 0.6;
    defaultOpts.widthBox = 0.3;
    defaultOpts.dodge = 0;
    defaultOpts.outliers = 0;
    defaultOpts.whiskers = 0;
    defaultOpts.notch = 0;
    defaultOpts.doJitter = 0;
    defaultOpts.ylim = [];
    defaultOpts.plotSize = [1 1 1600 900];
    defaultOpts.reorder = 0;
    defaultOpts.reorderX = 0;
    
    defaultOpts.reorderXbyY = 0;
    defaultOpts.npoints = 100;
    
    defaultOpts.title = [];
    defaultOpts.ylabel = '';
    defaultOpts.xlabel = '';
    defaultOpts.base_size = 12;
    
    defaultOpts.groupMode = 1;
    defaultOpts.dropLegend = 0;
    defaultOpts.cmap = [];
    defaultOpts.XTickLabelRotation = [];
    defaultOpts.equalizeY = 0;
    defaultOpts.orderFunc = @(x)nanmedian(x(:));
    
    defaultOpts.vNormalization = 'equalwdith';
    defaultOpts.vHistBin = 1;
    defaultOpts.bwfactor = 0.9;
    
    
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts)

    if isempty(opts.dodge)
        opts.dodge = 0;
    end
        
    zfig = [];
    if opts.newPlot == 1
        zfig = figure('Position',opts.plotSize);        
    end   
    
    [N,nC] = size(cX);
    
    % Alt match vector 
    if nargin<3 || isempty(cFact)
        cFact = [];                        
    else
        cFact = matchVector(cY,cFact);
    end
        
    cX = matchVector(cY,cX);            


        
    cX = cX(:);
    cY = cY(:);
    if ~isempty(cFact)
        cFact = cFact(:);
    end
    
    [cXlist,~,~,~,cXpos] = fastUnique(cX);
    if opts.reorderX == 1
        [cXlist,xReord] = sort(cXlist);
        
        cXpos = cXpos(xReord);        
    elseif isa(opts.reorderX,'containers.Map')
        [~,xReord] = sort(nanvalues(opts.reorderX,cXlist));
        cXlist = cXlist(xReord);
        cXpos = cXpos(xReord);
    elseif opts.reorderXbyY 
        cMVal = cellfun(@(x)opts.orderFunc(cY(x)),cXpos);
        
        [~,xReord] = sort(cMVal);
        
        cXpos = cXpos(xReord);        
        cXlist = cXpos(xReord);                
    end
    
    nC = length(cXlist);
    if opts.nSplit > 0
        opts.maxRow = min(opts.maxRow,ceil(nC/opts.nSplit));       
    end
    
    
    nR = ceil(nC/opts.maxRow);            
    zp = 1;    
    
    for i = 1:nR
        
        cP = zp:min(zp+opts.maxRow-1,nC);
        rowLen = length(cP);
        
        cSub = [ cXpos{cP} ];   
        
        if isempty(cFact)
            g(i,1) = gramm('x',cX(cSub),'y',cY(cSub));
        else
            g(i,1) = gramm('x',cX(cSub),'y',cY(cSub),'color',cFact(cSub));
        end
        
        if ~isempty(opts.widthV)

            g(i,1).stat_violin('width',opts.widthV,'dodge',opts.dodge,'npoints',opts.npoints)

        end
        if ~isempty(opts.widthBox)
            g(i,1).stat_boxplot('width',opts.widthBox,'dodge',opts.dodge,'outliers',opts.outliers,'whisk',opts.whiskers,'notch',opts.notch);
        end
        if opts.doJitter > 0
            g(i,1).geom_jitter('dodge',opts.dodge,'width',opts.doJitter);
        end


        if ~isempty(opts.ylim)
            g(i,1).axe_property('YLim',opts.ylim);
        end
        
        if opts.reorder == 0
            g(i,1).set_order_options('x',0);
        elseif opts.reorderXbyY == 1
            g(i,1).set_order_options('x',0);
%             if ~isempty(cFact)
%                 [zFactList,~,~,~,zFactOrder] = fastUnique(cFact(cSub));
%                 
%             else
%                 [zMlist,~,~,~,zMpos] = fastUnique();            
%                 cMedianY =   
%             emd
%             if isempty(cFact)
%                 
%                 [cXsub,~,~,~,cXsubPos] = fastUnique(x(cSub));
%                 
%                 cYsub = cY(cSub);
%                 cMVal = cellfun(@(x)median(cYsub(x)),cXsubPos);
%                 
%                 
%             else
%                 error('Not implemented');
%             end
%             
        end        
                  
        g(i,1).set_text_options('base_size',opts.base_size);
        if ~isempty(opts.title) && i == 1
            g(i,1).set_title(opts.title);            
        end

        
%         if ~isempty(opts.xlabel) && ~isempty(opts.xlabel)
%             
%         elseif ~isempty(opts.ylabel)
%             g(i,1).set_names('y',opts.ylabel);        
%         elseif ~isempty(opts.xlabel)
%             g(i,1).set_names('x','');
%         end
        g(i,1).set_names('x','','y',opts.ylabel);
        zp = zp + opts.maxRow;
    end
    
    if ~isempty(opts.xlabel)
        g(i,1).set_names('x',opts.xlabel,'y',opts.ylabel);
    end

    
    if opts.dropLegend == 1
        g.set_layout_options('legend',false);
    end
    if ~isempty(opts.cmap)
        g.set_color_options('map',opts.cmap,'n_color',size(opts.cmap,1));
    end    
           
    outG = g.draw();
    
    if opts.equalizeY
        cYlim = [];
        for zzc = 1:length(outG)
            cYlim(:,zzc) = outG(zzc).facet_axes_handles.YLim;
        end
        
        newYlim = [ min(cYlim(1,:)) max(cYlim(2,:)) ];
        
        for zzc = 1:length(outG)
            outG(zzc).facet_axes_handles.YLim = newYlim;
        end
        outG = outG.redraw();
    end
    
    
    if ~isempty(opts.XTickLabelRotation)
        for zzc = 1:length(outG)
            outG(zzc).facet_axes_handles.XTickLabelRotation = 45;
        end
        outG = outG.redraw();
    end
    
end

% function cY = matchVector(cX,cY,N,Nc)
% 
%     if ~all(size(cX) == size(cY))
%         if numel(cX) == numel(cY) && (min(size(cX)) == 1 || min(size(cY)) == 1)
%             warning('Number of elements matches but dimensions do not\n');
%         elseif size(cY,1) == N && size(cY,2) == 1
%             cY = repmat(cY,1,N);
%         elseif size(cY,2) == Nc && size(cY,1) == 1
%             cY = repmat(cY,N,1);
%         else
%             error('Size missmatch on %s',inputname(1));            
%         end            
%     end
% 
% end


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
