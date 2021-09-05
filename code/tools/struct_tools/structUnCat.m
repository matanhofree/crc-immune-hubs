function outCat = structToCat(strIn)

    if ~isstruct(strIn)
        error('Not a struct');
    end

    fList = fieldnames(strIn);
    
    outCat = strIn;
    for i = 1:length(fList)
        cf = fList{i};
        
        cV = strIn.(cf);
        
        if iscategorical(cV)
            outCat.(cf) = cellstr(cV);            
        elseif isstruct(cV)
            outCat.(cf) = structToCat(cV);
        end
    end

end