function [tselect,gout] = strgrep(strIn,regE,showSummary)
    
    if nargin < 3
        showSummary = 0;        
    end
    
    if iscategorical(strIn)
        strIn = cellstr(strIn);
    end
    if ~iscell(strIn)
        strIn = {strIn};
    end
       if isstring(regE)
        regE = cellstr(regE);
    end
    tt = regexp(strIn,regE);
        
    tselect = ~isemptycell(tt);
    gout = strIn(tselect);
    
    if showSummary && luniq(strIn) < 200
        tabFilter(strIn,tselect)
    end
end
    