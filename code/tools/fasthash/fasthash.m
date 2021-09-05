function [idxvector,hashKey,hN,binKey] = fasthash(hashKey,inVector)

    % hashKey = unique(hashKey);
    if luniq(hashKey) ~= length(hashKey)
        warning('Key is not unique');
        hashKey = fastUnique(hashKey);
    end
    
    if strcmp(class(hashKey),class(inVector)) ~= 1
        error('Input missmatch');
    end      
    
    if iscell(hashKey)    
        hashKey(isemptycell(hashKey)) = [];    
        zdrop = ~isemptycell(inVector);
    
        idxvector = nan(length(inVector),1);
        idxvector(zdrop) = (fasthashStrPtr(hashKey,inVector(zdrop)));
    else
        %hashKey(isempty(hashKey)) = [];    
        % zdrop = ~isempty(inVector);
    
        % idxvector = nan(length(inVector),1);
        idxvector = fasthashInt(hashKey,inVector);
    end
        
    if nargout > 2
        [hN,binKey] = histcounts(idxvector,1:length(hashKey)+1);
        hN = hN(:);
    end
end