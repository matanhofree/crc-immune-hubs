function auc = cAUC(predVal,rlabel)      
% Calculate AUC
% AUC: a Statistically Consistent and more Discriminating Measure than Accuracy, by Charles X. Ling, Jin Huang and Harry Zhang
% 
      zdrop = isnan(predVal(:)) | isnan(rlabel(:));
      predVal(zdrop) = [];
      rlabel(zdrop) = [];

      nTarget     = sum(double(rlabel == 1));
      nBackground = sum(double(rlabel == 0));

      % Rank data
      R = tiedrank(predVal);  % 'tiedrank' from Statistics Toolbox

      % Calculate AUC
      auc = (sum(R(rlabel == 1)) - (nTarget^2 + nTarget)/2) / (nTarget * nBackground);
