function [tselect,gout] = strgrepE(strIn,regE)
    
    if ~iscell(strIn)
        strIn = {strIn};
    end
    
    if numel(regE) > 1 && iscellstr(regE)
        regE = strjoin(regE,'|');
    end
    tt = regexp(strIn,regE);
        
    tselect = ~isemptycell(tt);
    gout = strIn(tselect);
    
end
    