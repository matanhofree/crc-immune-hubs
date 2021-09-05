function [zfigOut,outTables,outMat] = ccNMF_plot_coAct(outData,outPath,inOpts)

    
    defaultOpts.sigThr = [];
    defaultOpts.sigThrTxt = [];
    defaultOpts.cropCellName = 1;
    defaultOpts.annotCell = 1;
    defaultOpts.rpower = 1;
    
    defaultOpts.saveXlsx = 1;
    
    defaultOpts.outTables = 4;
    defaultOpts.suffixPlot = { '-dpng' };
    defaultOpts.clobber = 0;
    defaultOpts.doFisher = 1;
    
    defaultOpts.doSub = 1;

    defaultOpts.pdist = 'corr';
    defaultOpts.linkage = 'average';
    

    defaultOpts.coAct.showTxt = 0;
    defaultOpts.coAct.maxLabelX = 150;
    defaultOpts.coAct.maxLabelY = 150;
    % defaultOpts.coAct.colorHard = [ -1 1 ];
    defaultOpts.doTrimNames = 1;
    defaultOpts.clMap = [];
    
    defaultOpts.sigRefFdr = [];

    


    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts
    
    outTablesSuffix = { 'meanR' 'Radj' 'FDR' 'Pval' };
    doSave = 0;
    
    if nargin > 1 && ~isempty(outPath)
        doSave = 1;
        [outPathBase,outName] = fileparts(outPath);
        
    end
        
    
    zfigOut = [];
    outTables = [];
    outMat = [];
    
    optsCoAct = opts.coAct;
        
    outC = outData.outMS;
    annotName = outC.hMatNames;
    
%     trimNode = @(x)regexprep(x,'.*_zz_','');
%     if opts.doTrimNames
%        annotName = trimNode(annotName);
%     end

    
    if all(strgrep(annotName,'^p'))
        cellType = regexprep(annotName,'p([^0-9]*).*','$1');
    else
        cellType = regexprep(annotName,'(.*)_zz_.*','$1');
    end

    if opts.annotCell
        optsCoAct.annotX.sampleID = annotName;
        optsCoAct.annotX.cellType = cellType;
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
       %  adjR = max(min(adjR,1),-1);
        adjR(eye(size(adjR,1))>0) = 0;

    end
    
    zc = sum(adjR > 1);
    if zc > 0
        fprintf('Warning: %d values > 1 in adjR',zc);
    end
    zc = sum(adjR < -1);
    if zc > 0
        fprintf('Warning: %d values < -1 in adjR',zc);
    end
   
    
    cDistR = adjR + eye(size(outC.corrP,1));
    % cDistR = max(min(cDistR,1),-1);

%     cDistR = cDistR - min(cDistR(:));
%     cDistR = cDistR - min(cDistR(:));
    
    cDistR = squareform(1 - cDistR).^opts.rpower;
    optsCoAct.doLeafOptimalOrderX_data.pdist = cDistR;
    optsCoAct.doLeafOptimalOrderY_data.pdist = cDistR;
    
    %%
%     cPval = min(min(outC.corrCountP,outC.corrCountN)*2,1);
%     cFdr = cPval;
    nPv = size(outC.corrCountP,1);
    
    cFdr = nan(nPv);
    cPval =  min(min(outC.corrCountP,outC.corrCountN)*2,1);
    
    zSelU = triu(true(nPv),1);
    
    cFdrPre = [ outC.corrCountP(zSelU) outC.corrCountN(zSelU) ];
    cFdrPre(:) = mafdr(cFdrPre(:),'BHfdr',1);
    %%

%     if ~issymmetric(cPval)
%         zz = tril(cPval,-1) - triu(cPval,1)';
%         fprintf('Warning: not symmetric. Diff %d, max diff',sum(zz(:)),max(zz(:)));
%     end    %%
    cFdr(zSelU) = min(cFdrPre,[],2);
    tmpT = cFdr';
    cFdr(zSelU') = tmpT(zSelU');
    % assert(issymmetric(cFdr));
    
    %%
    adjRplot = adjR;
    if ~isempty(opts.sigThr)
        cSelSig = cFdr < opts.sigThr;
        
        adjRplot(~cSelSig) = 0;
        
        if sum(cSelSig(:)) < 5
            fprintf('Insufficient output to plot');
            return;
        end
    end
    %%

    sigTxt = [];
    if ~isempty(opts.sigThrTxt)
%%
        sigTxt = cell(size(adjR,1));
        sigTxt(:) = {''};
        optsCoAct.showTxt = 2;
%%
        if ~isempty(opts.sigRefFdr)
            %%
            fprintf('Using external reference FDR\n');
            
            zSelU = triu(true(nPv),1);
            cFdrExt = nan(nPv);
            
            [~,zia,zib] = intersect(annotName,opts.sigRefFdr.id);
            
            extR = nan(nPv);
            extFDR = nan(nPv);                        
            
            extR(zia,zia) = opts.sigRefFdr.adjR(zib,zib);
            extFDR(zia,zia) = opts.sigRefFdr.cFdr(zib,zib);
                        
            cSelSubU = zSelU(:) & extR(:) > 0;% & extFDR(:) < 0.1;
            cSelSubD = zSelU(:) & extR(:) < 0;% & extFDR(:) < 0.1;

            % cSelSubUref = extR(:) > 0 & extFDR(:) < 0.1;
            % cSelSubDref = extR(:) < 0 & extFDR(:) < 0.1;
            cSelSubRef = extFDR(:) < 0.1;
                        
            cPvalRef = [ outC.corrCountP(cSelSubU); outC.corrCountN(cSelSubD) ];
            cFdrRef = mafdr(cPvalRef,'BH',1);
            
            cFdrExt(cSelSubU) = cFdrRef(1:sum(cSelSubU));
            cFdrExt(cSelSubD) = cFdrRef(sum(cSelSubU)+1:end);
                    
            
            tmpT = cFdrExt';
            cFdrExt(zSelU') = tmpT(zSelU');
            
            sigTxt(cSelSubRef) = {'o'};
            sigTxt(cFdrExt < opts.sigThrTxt) = {'*'};
            
        else
            cSelSig = cFdr < opts.sigThrTxt;

            sigTxt(cSelSig) = {'*'};
        end
    end
    
    if ~isempty(opts.clMap)
%         optsCoAct.moduleSortX = 1;
%         optsCoAct.moduleSortY = 1;
        optsCoAct.doLeafOptimalOrderX = 0;
        optsCoAct.doLeafOptimalOrderY = 0;
        optsCoAct.doSortY = 0;
        optsCoAct.doSortX = 0;

        
        modMap = containers.Map(opts.clMap.nodeIdx,opts.clMap.moduleIdx);
        
        clV = nanvalues(modMap,annotName);
        if any(isemptycell(clV))
            cNan = isemptycell(clV)
            cNanName = mergeStringPair('zz',1:sum(cNan));
            clV(cNan) = cNanName;
        end
        
        if ~isfield(opts.clMap,'leafOrder')
            leafOrder = constrainedLeafOrder(adjRplot,clV,2,opts.pdist,opts.linkage,1,opts);        
        else
            fprintf('Using input order');
            [~,zia,zib] = intersect(annotName,opts.clMap.leafOrder);
            
            leafOrder = zia(argsort(zib));
            
        end
            
        optsCoAct.externXOrder = leafOrder;
        optsCoAct.externYOrder = leafOrder;
        
        optsCoAct.annotX.modules = clV;

    end
    

    [zfigOut,outMat,outMat.zOrdTX,outMat.zOrdTY] = plot_heatmap_annot(adjRplot,annotName,annotName,[],sigTxt,optsCoAct);
    
    if doSave   
        if ~isempty(opts.sigThr)
            cOutF = [outName '_sigOnly_ccHeatMap' ];
        else
            cOutF = [outName '_ccHeatMap' ];
        end
        
        print_plot(zfigOut,cOutF,outPathBase,opts.suffixPlot,opts.clobber);
    end
    
    
    outTables = [];    
    cTableName = matlab.lang.makeValidName(annotName(outMat.xOrder));
    
    if opts.outTables >= 1
        zCC = outC.corrP;
        
        outTables{1} = [ cell2table(annotName(outMat.yOrder),'VariableNames',{ 'id' }) ... 
                      array2table(zCC(outMat.yOrder,outMat.xOrder),'VariableNames',cTableName) ];
   
    end
    if opts.outTables >= 2
        zCC = adjR;
        
        outTables{2} = [ cell2table(annotName(outMat.yOrder),'VariableNames',{ 'id' }) ... 
                      array2table(zCC(outMat.yOrder,outMat.xOrder),'VariableNames',cTableName) ];
               
    end
    if opts.outTables >= 3
        zCC = cFdr;
        
        outTables{3} = [ cell2table(annotName(outMat.yOrder),'VariableNames',{ 'id' }) ... 
                      array2table(zCC(outMat.yOrder,outMat.xOrder),'VariableNames',cTableName) ];
                             
    end
    if opts.outTables >= 4
        zCC = cPval;
        
        outTables{4} = [ cell2table(annotName(outMat.yOrder),'VariableNames',{ 'id' }) ... 
                      array2table(zCC(outMat.yOrder,outMat.xOrder),'VariableNames',cTableName) ];
         
    end
    
    if doSave
        for ci = 1:length(outTables)
            if isempty(outTables{ci})
                continue;
            end
                    
            if opts.saveXlsx   
                writetable(outTables{ci},[ outPath '.xlsx'],'sheet',outTablesSuffix{ci});
            else
                writeTableFile(outTables{ci},sprintf('%s_%s',outPath,outTablesSuffix{ci}),opts.clobber);
            end   
        end
    end
    
end