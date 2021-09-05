function [cM,outSetNames,mQ,opts] = graph_cluster_hc(G,fName,inOpts)
    
    defaultOpts.tempDir = 'temp/siGraphCl_temp'    
    defaultOpts.jarPath = '/Users/mhofree/projects/cancer_SC/temp/signed-community-detection/jar/simap-1.0.0-final.jar'    
    
    defaultOpts.jarParam = ' -i  0.00001 0.2 -a 0.0005 --tau 0.15';   
    
    defaultOpts.makeTemp = 1;
        
    defaultOpts.doFisher = 1;
    defaultOpts.upperOnly = 1;
    
    defaultOpts.posOnly = 0;

    defaultOpts.clearTemp = 0;
    
    defaultOpts.doPrint = 0;

    defaultOpts.minHC = 3;
    defaultOpts.minSubQ = 0;
    defaultOpts.optimizeSubQ = 0;

    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts);
    
    cM = []
    outSetNames = [];
    mQ = -inf;
    
    nN = size(G.Nodes,1);
    if nN < opts.minHC      
        return;
    end
    
    if ~exist('fName','var')
        fName = [];
    end

        
    [cM,outSetNames,mQ,oRes,rOpts] = graph_cluster_signed(G,fName,opts);
    
    if isnumeric(cM.moduleIdx)
        cM.moduleIdx = num2cellstr(cM.moduleIdx);
    end

    if luniq(cM.moduleIdx) == 1
        return; 
    end
    cPre = cM;
    [cModList,~,~,cModCount,cModPos] = fastUnique(cM.moduleIdx);
    
    for zi = 1:length(cModList)
        selSubIdx = trueV(cModPos{zi},nN);
        if cModCount(zi) > opts.minHC
            
            subName = sprintf('%s/%s_sub_%s',rOpts.fDir,rOpts.fName,cModPos{zi});
                
                
            subG = subgraph(G,selSubIdx);
            
            [subM,~,subQ] = graph_cluster_hc(subG,subName,opts);
                        
            if luniq(subM.moduleIdx) > 1 && subQ(1) > opts.minSubQ
                cM.moduleIdx(selSubIdx) = mergeStringPair(cM.moduleIdx(selSubIdx),subM.moduleIdx)                 
            end            
        end
    end
    
    if ~isequal(cPre,cM) && ~isempty(opts.optimizeSubQ)
    
        [cModList,~,~,cModCount,cModPos] = fastUnique(cPre.moduleIdx);
        cAdj = adjacency(G,'weighted');
        cQ = calcQ(cAdj,cPre.moduleIdx);
        finalIdx = cPre.moduleIdx;
        
        for zi = 1:length(cModList)
                        
            selSubIdx = trueV(cModPos{zi},nN);
            submoduleIdx = cPre.moduleIdx;
            submoduleIdx(selSubIdx) = cM.moduleIdx(selSubIdx);
            % submoduleIdx(selSubIdx) = regexprep(submoduleIdx(selSubIdx),'_[^_]*$','');
           
            subQ = calcQ(cAdj,submoduleIdx);
            cDiff = (cQ - subQ);
            if cDiff <=  opts.optimizeSubQ
                fprintf('Dropped %s -- %g\n',cModList{zi},cDiff);
            else
                finalIdx(selSubIdx) = submoduleIdx(selSubIdx);
                fprintf('Keeping %s -- %g\n',cModList{zi},cDiff);
            end
        end
        
        cM.moduleIdx = finalIdx;
        [cModList,~,~,cModCount,cModPos] = fastUnique(cM.moduleIdx);
        nodeNames = G.Nodes.Names;
        for i = 1:length(cModList)
        
            outSetNames{i} = nodeNames(cModPos{i});
        end
    end
    
    if ~isempty(outSetNames)
        outSetNames = cell2struct(outSetNames',mergeStringPair('%s%s','c',cModList));
    end

end


function [subQ,mQ] = calcQ(cAdj,cMod)
    cAdjP = cAdj;
    cAdjP(cAdj < 0) = 0;
    cAdjN = abs(cAdj);
    cAdjN(cAdj > 0) = 0;
    
    mQ(2:3) = [ weightedModularityQ(cAdjP,cMod) -weightedModularityQ(cAdjN,cMod) ];
    mQ(1) = nansum(mQ);
    subQ = mQ(1); 
    
    
end