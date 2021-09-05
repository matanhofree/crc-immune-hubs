function [isemptyarray] = isemptycell(cellarray)
    
    isemptyarray = cellfun(@(x)isempty(x),cellarray);
end