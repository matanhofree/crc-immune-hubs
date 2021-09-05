function [outTopIdx,outTopGenes,outTopMerged,outTopOrder,outTopOrderGene] = ccNMFextractTopGenes(W,geneID,dynamicN,hardN)

    Wmi = weightedNNcontrast(W);    
    
    [D,K] = size(Wmi);
        
    [cVal,cidxWord] = sort(Wmi,'descend','MissingPlacement','last');
        
    % 
    if isempty(dynamicN)
        cValThr = cVal(hardN,:);

    else
        fprintf('Using dynamic elbow based threshold rule.\n');
        cValThr = findElbow(cVal(1:dynamicN,:)) ;
        cValThr = max(cValThr,cVal(hardN,:));
    end
%     else 
%         fprintf('Using hard rank threshold rule.\n');
%         
%         cValThr = cVal(opts.geneTableTopN,:);
%     end
%%
%    cValThr(cValThr == 0) = nan;
    outTopIdx = Wmi > cValThr;
    
    outTopIdx(isnan(Wmi)) = 0;    
    
    fprintf('Found genes:');  
    cidxGeneSelectTotal = sum(outTopIdx);
    disp(cidxGeneSelectTotal);
  %%          
    if ~isempty(geneID)
        outTopGenes = cell(max(cidxGeneSelectTotal),K);
        outTopOrder = cell(K,1);
        outTopOrderGene = cell(K,1);
        
        outTopGenes(:) = { '' };
    
        for i = 1:K
            cSel = 1:cidxGeneSelectTotal(i);
            cG = geneID(cidxWord(cSel,i));
            outTopGenes(cSel,i) = cG;
            outTopOrderGene{i} = cG;
            
            outTopOrder{i} = cidxWord(cSel,i);
        end
    else
        outTopGenes = [];
    end
    
    outTopMerged = any(outTopIdx,2);
    fprintf('Identified a total of %d genes\n',sum(outTopMerged));
    
end
