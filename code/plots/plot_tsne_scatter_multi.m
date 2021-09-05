function zfig = plot_tsne_scatter_multi(ydata,mVar,colorSet,sliceVal,inOpts)

    defaultOpts.smpNames = [];    
        
    defaultOpts.nCol = 4;
    defaultOpts.nRow = 3;       
    defaultOpts.plotSize = [ 10 10 2010 1510 ];
    defaultOpts.tightAxis = 1;
    defaultOpts.isCategory = 0;
    defaultOpts.applyZscore = 0;
    defaultOpts.qTrimLim = 0.01;
    defaultOpts.qTrimLimNNZ = 0.99;
    defaultOpts.cLimMinMax = 1;
    defaultOpts.cLimHard = [];
    
    defaultOpts.titleText = 'C %d';
    defaultOpts.tiledLayout = 0;
    

    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    
    if nargin < 3
        colorSet = [];
    end
    
    if nargin < 4
        sliceVal = [];
    end

    nV = size(ydata,1);
    
    if isstruct(mVar)
        cf = fieldnames(mVar);
        mVar = struct2cell(mVar);
        isSmp = strgrep(cf,'sampleID');
        if any(isSmp)
            mVar(isSmp) = [];
            cf(isSmp) = [];
        end 
        
        opts.titleText = cf;
    end
    if iscell(mVar)
        mLen = cellfun(@(x)length(x),mVar);
        if ~all(mLen == nV)
            error('Dim missmatch');
        end
    else 
        if size(mVar,2) == nV
            mVar = num2cell(mVar,2);        
        elseif size(mVar,1) == nV
            mVar = num2cell(mVar,1);        
        else
            error('Dim missmatch');
        end        
    end
    
    zfig = [];
    
    plotPerFig = opts.nCol*opts.nRow;
    
    nV = length(mVar);
    zn = 1;
    zf = 1;
    for zi = 1:nV
        cStat = mVar{zi};
        if ~iscell(cStat)
            if all(isnan(cStat))
                if ~isempty(opts.titleText)
                    if iscell(opts.titleText)
                        cSt = sprintf(opts.titleText{zi});
                    else
                        cSt = sprintf(opts.titleText,zi);
                    end
                    fprintf('Skipping %s\n',cSt);            
                else
                    fprintf('Skipping %d\n',zi);
                end            
                continue;
            end            
        end        
        if zn > plotPerFig || zf == 1
           zfig{zf} = figure('Position',opts.plotSize,'color','w');
           zf = zf+1;
           zn = 1;
           
           if opts.tiledLayout
               tiledlayout(opts.nRow,opts.nCol,'Padding', 'compact', 'TileSpacing', 'compact'); 
           end
        end
        
        if opts.tiledLayout
           	nexttile();            
        else
            zaxes = subplot(opts.nRow,opts.nCol,zn);             
        end
        zn = zn+1;

        zopts = opts;
        zopts.newPlot = 0;
        zopts.cmapQuant = 0;
        
        if  ~iscell(cStat) && opts.applyZscore
            cStat = zscore(cStat);
            zopts.yAxis_unit = 'z-score';
        end 
        
        if ~iscell(cStat)
            if isempty(opts.cLimHard)
                if min(cStat) < 0 && ~isempty(opts.qTrimLim)
                    zopts.cmapLimits = quantile(cStat,[opts.qTrimLim 1-opts.qTrimLim]);
                elseif ~isempty(opts.qTrimLimNNZ);
                    zopts.cmapLimits = [ 0 quantile(nonzeros(cStat),opts.qTrimLimNNZ) ];              
                elseif opts.cLimMinMax 
                    zopts.cmapLimits = [ min(cStat) max(cStat) ];
                end            
            else
               zopts.cmapLimits = zopts.cLimHard; 
            end
        
        else
            zopts.isCategory = 1;
        end
        
        plot_tsne_scatter(ydata,cStat,colorSet,sliceVal,zopts);
        
        if opts.tightAxis
            axis tight 
            axis off 
        end
        
        if ~isempty(opts.titleText)
            if iscell(opts.titleText)
                title(opts.titleText{zi},'Interpreter','none');
            else
                title(sprintf(opts.titleText,zi));
            end
        end           
    end   
end