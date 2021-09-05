function stOut = structRowsToColumns(stIn)

    zn = fieldnames(stIn);
    
    for i = 1:length(zn)
        cf = zn{i};
        [xi,xj] = size(stIn.(cf));
        if (xi == 1)
            stOut.(cf) = stIn.(cf)';
        else
            stOut.(cf) = stIn.(cf);
        end
    end
    