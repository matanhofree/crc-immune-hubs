function [figOut,outMat,zOrdTX,zOrdTY] = plot_crosstab_heatmap_alt(clustMain,clustStack,cmap,inputTable,inOpts)

    defaultOpts.axisObj =[];
    defaultOpts.clustTxt = 1;
    
    defaultOpts.XTickLabelRotation = 60;
    
    defaultOpts.colormap = flipud(cptcmap('ylwhbl_sym'));
    defaultOpts.maxColorLim = 0.95;
    defaultOpts.fontSize = 18;
    defaultOpts.axisFontSize = 18;

    defaultOpts.pseudoCount = 5;
    defaultOpts.colorType = 'pearson';
    
    defaultOpts.annotX = [];
    defaultOpts.extPR = [];
    
    defaultOpts.gridLines = 'none';
            
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts);
    zfig = [];
    
    if ~exist('cmap','var') || isempty(cmap)
        cmap = opts.colormap;
    end
   
    if exist('inputTable','var') && ~isempty(inputTable)
        crossTable = table2array(inputTable);
        xName = inputTable.Properties.VariableNames;
        yName = inputTable.Properties.RowNames;
    else
        [crossTable,~,~,headerC] = crosstab(clustMain,clustStack);
    
        yName = headerC(:,1);
        yName(isemptycell(yName)) = [];

        xName = headerC(:,2);
        xName(isemptycell(xName)) = [];
    end  


    % pTable = crossTable + opts.pseudoCount; 
    if isempty(opts.extPR)
        switch opts.colorType
            case 'pearson'
                fprintf('Using pearson residual transformation for table');
                pTable = crossTable + opts.pseudoCount; 
                cTotal = sum(pTable(:));
                pTableN = pTable./cTotal;    
                cExp = (sum(pTableN,2)*sum(pTableN))*cTotal;
                cPR = (pTable-cExp)./sqrt(cExp);
            % case 'pearsonNormalCol'
            %     fprintf('Using pearson residual transformation for table');

            %     pTableN = crossTable + opts.pseudoCount; 

            case 'logOdds'
                fprintf('Using log odds');
                pTable = crossTable + opts.pseudoCount; 
                cTotal = sum(pTable(:));
                pTableN = pTable./cTotal;    
                cExp = (sum(pTableN,2)*sum(pTableN))*cTotal;
                
                cPR = log2(pTable./cExp);
                
            case 'rowFreq'
                fprintf('Using rowFreq');
                pTable = crossTable + opts.pseudoCount; 
                cPR = pTable./sum(pTable,2);
                
            case 'colFreq'
                fprintf('Using colFreq');
                pTable = crossTable + opts.pseudoCount; 
            
                cPR = pTable./sum(pTable,1);
            otherwise 
                error('Pick coloring by - pearson or logOdds');                
        end    
    else  
        fprintf('Using external PR');
        cPR = opts.extPR;

        if any(size(cPR) ~= size(crossTable))
            error('Missmatch between input PR table and crossTable found.');
        end
    end

    if ~isempty(opts.annotX)
        %%
        flist = setdiff(fieldnames(opts.annotX),'sampleID');
        [xNameList,~,~,~,xNamePos] = fastUnique(opts.annotX.sampleID);
         
        reAnnotX.sampleID = xName;
        for i = 1:length(flist)
            cf = flist{i};
            cV = opts.annotX.(cf);
            % if length(cV) ~= length(clustStack)
            %     error('Annotation size missmatch');
            % end
            cAnnotMap = containers.Map(xNameList,cellfun(@(x)strjoin(unique(cV(x)),';'),xNamePos,'unif',0));
            
            outV = nanvalues(cAnnotMap,xName);
            outV(isemptycell(outV)) = { 'NA' };
            reAnnotX.(cf) = outV; 
        end   
        opts.annotX = reAnnotX;
    end
    %%
    [figOut,outMat,zOrdTX,zOrdTY] = plot_heatmap_annot(cPR,yName,xName,cmap,crossTable,opts);       
    
end
