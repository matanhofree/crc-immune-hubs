function [ meanExp, fracExp, batchNames, outAUC, outPRAUC, countB ] = summarize_subset_value(aData,aBatch,inOpts)
    
    defaultOpts.minA = [];

    defaultOpts.applyLogA = 0;
    
    defaultOpts.doLog = 0;
    
    defaultOpts.dimR = 2;
    defaultOpts.aggrFunc = @(X,dim)nanmean(X,dim);
    defaultOpts.nFrac = @(X,dim)nansum(X>0,dim)/size(X,dim);
    
    defaultOpts.doFigures = 1;
    defaultOpts.doSort = 1;
    defaultOpts.doNaNzeros = 0;
    
    defaultOpts.doAUC = 1;
    defaultOpts.doPRAUC = 0;    
    defaultOpts.pb = 1;
        
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    % disp(opts);   
    
    if (nargout < 4 && nargout > 1)
        opts.doAUC = 0;        
    end
    
    if (nargout < 5 && nargout > 1)
        opts.doPRAUC = 0;        
    end        
        
    aggrFunc = @(X)opts.aggrFunc(X,opts.dimR);
    nFrac = [];
    if ~isempty(opts.nFrac)
        nFrac = @(X)opts.nFrac(X,opts.dimR);
    end
      
    [batchNames,~,~,countB,cPos] = fastUnique(aBatch);  
    if opts.doSort
        [~,bIdx] = sort(batchNames);    
        batchNames = batchNames(bIdx);
        countB = countB(bIdx);
        cPos = cPos(bIdx);
    end
    
    nBatchA = length(batchNames);
    
    [nD,nN] = size(aData);
    
    outAUC = [];
    outPRAUC = [];
    fracExp = [];
    if opts.dimR == 2
        meanExp = nan(nD,nBatchA);
        if ~isempty(nFrac)
            fracExp = nan(nD,nBatchA);
        end
        if opts.doAUC
            outAUC = nan(nD,nBatchA);
        end           
        if opts.doPRAUC
            outPRAUC = nan(nD,nBatchA);
        end        
    elseif opts.dimR == 1
        meanExp = nan(nN,nBatchA);
        if ~isempty(nFrac)
            fracExp = nan(nN,nBatchA);
        end
        if opts.doAUC
            error('Not yet implemented');
            outAUC = nan(nN,nBatchA);
        end           
        if opts.doPRAUC
            error('Not yet implemented');
            outPRAUC = nan(nN,nBatchA);
        end
    end
        
    if opts.pb && ismac()
        progressbar();
    end
    
    for zi = 1:nBatchA
%         if iscell(aBatch)            
%             zIdx = strcmp(aBatch,batchNames{zi});        
%         else
%             zIdx = (aBatch == batchNames(zi));        
%         end
        zIdx = trueV(cPos{zi},nN);
            
        aDataSub = aData(:,zIdx);
        
        if opts.doNaNzeros == 1
            fprintf('Zeros to nans\n');
            aDataSub = full(aDataSub);
            
            aDataSub(aDataSub == 0) = nan;
        end
        
        zEst = aggrFunc(aDataSub);
        if ~isempty(nFrac)
            zFrac = nFrac(aDataSub);
        end
        
        if ~isempty(opts.minA)
            zEst(zEst<opts.minA) = nan;
        end
        
        if opts.applyLogA == 1
            zEst = log1p(zEst);
        end
        meanExp(:,zi) = zEst;
        if ~isempty(nFrac)
            fracExp(:,zi) = zFrac;
        end
        
        if opts.doAUC
            outAUC(:,zi) = arrayfun(@(xi)cAUC(aData(xi,:),zIdx),1:nD);
        end
        
        if opts.doPRAUC
            for zj = 1:nD
                [~,~,~,outPRAUC(zj,zi)] = perfcurve(zIdx,aData(zj,:),1,'Xcrit','prec','YCrit','reca');
            end
        end
        if opts.pb && ismac()
            progressbar(zi/nBatchA);
        end
    end  
    
    if nargout == 1
        outD = structpack([],meanExp,fracExp,batchNames,outAUC,outPRAUC,countB,cPos);
        meanExp = outD;
    end
    
end