function wq = weightedModularityQ(aMat,modCl)

    wR = sum(aMat,2);
    wC = sum(aMat);
    wS = 2*sum(wR);
    
    %%
    gB = aMat - ((wR*wC)/wS);
    
    %%
    [clB,~,~,~,clPos] = fastUnique(modCl);

    gS = zeros(size(aMat,1),length(clB));
    for i = 1:length(clB)
        gS(clPos{i},i) = 1;
    end
    
    %%
    clQ = gS'*gB*gS;
    wq = trace(clQ)/wS;

%% 
%     
%     zG = gS(:,1) - gS(:,2);
% %     
% %     
%      zz = (zG'*gB*zG)/(wS)
% %     
end
