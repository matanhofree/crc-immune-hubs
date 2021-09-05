function [adjR,cPval,cFdr,rawR] = extractSquarePvalFdr(outC,inOpts,extR)

    defaultOpts.doSub = 1;
    defaultOpts.doFisher = 1;
    defaultOpts.useRsign = 0;
    
    

    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts
   
    rawR = outC.corrP;
    
    if nargin < 3
        extR = [];
    end
    
    if opts.doSub == 1
        if opts.doFisher == 1
            adjR = tanh(atanh(outC.corrP) - atanh(outC.permMean));
            adjR(eye(size(adjR,1))>0) = 0;
        else 
            adjR = outC.corrP - outC.permMean;
            adjR = max(min(adjR,1),-1);
        end
    else
        adjR = outC.corrP;
        adjR(eye(size(adjR,1))>0) = 0;
    end
   
    nPv = size(outC.corrCountP,1);
    if opts.useRsign == 0
        fprintf('Simple p-value\n');
        
        cFdr = nan(nPv);
        cPval =  min(min(outC.corrCountP,outC.corrCountN)*2,1);

        zSelU = triu(true(nPv),1);

        cFdrPre = [ outC.corrCountP(zSelU) outC.corrCountN(zSelU) ];
        cFdrPre(:) = mafdr(cFdrPre(:),'BHfdr',1);
        %%


        cFdr(zSelU) = min(cFdrPre,[],2);
        tmpT = cFdr';
        cFdr(zSelU') = tmpT(zSelU');
        % assert(issymmetric(cFdr));
    else
        if ~isempty(extR)
            fprintf('Using external R to pick direction\n');
            extR = rawR;
        else
            fprintf('Using rawR to pick direction\n');
            extR = rawR;
        end
            
        cPval = nan(nPv);               
        cSel = extR > 0 & adjR > 0;
        cPval(cSel) = outC.corrCountP(cSel);
        
        cSel = extR < 0 & adjR < 0;
        cPval(cSel) = outC.corrCountN(cSel);
        
        cFdr = cPval;
        zSelU = triu(true(nPv),1);

        zSelV = zSelU(:) & ~isnan(cFdr(:));
        cFdr(zSelV) = mafdr(cFdr(zSelV),'BHfdr',1);
        
        tmpT = cFdr';
        cFdr(zSelU') = tmpT(zSelU');
    end
    % zz = [ adjR(:) cPval(:) outC.corrCountP(:) outC.corrCountN(:) rawR(:)]
end