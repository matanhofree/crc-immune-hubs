function outT = cellToColTable(inCell,inF)
    
    if nargin < 2
        inF = [];
    end

    if isstruct(inCell)
        inF = fieldnames(inCell);
        inCell = struct2cell(inCell);
    end
        

    cD = length(inCell);
    cN = cellfun(@(x)numel(x),inCell);
    
    outT = cell(max(cN),cD);
    outT(:) = { '' };
    
    for i = 1:cD
        outT(1:cN(i),i) = inCell{i};
    end
    
    if isempty(inF)
        outT = cell2table(outT);
    else
        outT = cell2table(outT,'variableNames',inF);
    end
end