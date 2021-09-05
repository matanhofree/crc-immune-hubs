function outCat = structToCat(strIn,maxT)

    if ~isstruct(strIn)
        error('Not a struct');
    end

    if nargin < 2
        maxT = 100;
    end

    fList = fieldnames(strIn);
    
    outCat = strIn;
    for i = 1:length(fList)
        cf = fList{i};
        
        cV = strIn.(cf);
        
        if iscell(cV) && luniq(cV) < maxT
            outCat.(cf) = categorical(cV);            
        end
    end

end