function [outT,chi,p,zRowH,zColH] = crossTabTable(varX,varY,normD,sortIdx,pseudoB)

    if nargin < 3
        normD = 0;
    end

    if nargin < 4
        sortIdx = 0;
    end

    if nargin < 5
        pseudoB = 0;
    end


    
    [zTab,chi,p,zTabH] = crosstab(varX,varY);
    
    zTab = pseudoB + zTab;
    
    switch normD
        case 1 
            zTab = zTab./sum(zTab,1);
        case 2
            zTab = zTab./sum(zTab,2);
        case 3
            zTab = zTab/sum(zTab(:));
        case 4
            zN = sum(zTab(:));
            zTabE = zTab/zN;
            zExp = (sum(zTabE).*sum(zTabE,2))*zN;            
            zTab = (zTab - zExp)./sqrt(zExp);
        case 5
            zN = sum(zTab(:));
            zTabE = zTab/zN;
            zExp = (sum(zTabE).*sum(zTabE,2))*zN;            
            zTab = zTab.*log(zTab./zExp);
        case 6
            zTab = zTab + 1;
            zN = sum(zTab(:));
            zTabE = zTab/zN;
            zExp = (sum(zTabE).*sum(zTabE,2))*zN;            
            zTab = log2(zTab./zExp);
    end
    
    zRowH = zTabH(:,1);
    zColH = zTabH(:,2);
    zRowH(isemptycell(zRowH)) = [];
    zColH(isemptycell(zColH)) = [];
    
    zRowH = matlab.lang.makeValidName(zRowH);
    zColH = matlab.lang.makeValidName(zColH);
    
        
    switch sortIdx
        case 1 
            [zRowH,zir] = sort(zRowH);
            zTab = zTab(zir,:);
        case 2
            [zColH,zic] = sort(zColH);
            zTab = zTab(:,zic);
        case 3
            [zRowH,zir] = sort(zRowH);
            zTab = zTab(zir,:);
            [zColH,zic] = sort(zColH);
            zTab = zTab(:,zic);
    end
    
    outT = array2table(zTab,'VariableNames',zColH,'RowNames',zRowH);
    
end