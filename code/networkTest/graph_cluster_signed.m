function [cM,outSetNames,mQ,oRes,opts] = graph_cluster_signed(G,fName,inOpts)
    
    defaultOpts.tempDir = 'temp/siGraphCl_temp'    
    defaultOpts.jarPath = '/Users/mhofree/projects/cancer_SC/temp/signed-community-detection/jar/simap-1.0.0-final.jar'    
    
    defaultOpts.jarParam = ' -i  0.001 0.2 --tau 0.1 ';    
    defaultOpts.makeTemp = 1;
        
    defaultOpts.doFisher = 1;
    defaultOpts.upperOnly = 1;
    
    defaultOpts.posOnly = 0;

    defaultOpts.clearTemp = 0;
    
    defaultOpts.doPrint = 1;

    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts);
        
    fDir = [];
    fName = [];
    if exist('fName','var') && ~isempty(fName)
        [fDir,fName] = fileparts(fName);        
    end

    if isempty(fDir)
        fDir = opts.tempDir;
    end 
    if isempty(fName)
        [~,fName] = fileparts(tempname);
    end
    
    if opts.makeTemp == 1
        mkdir(fDir);
    end
    
    opts.fDir = fDir;
    opts.fName = fName;
    
    cAdj = adjacency(G,'weighted');
    
    [zi,zj,zv] = find(cAdj);
    
    if opts.upperOnly 
        if ~issymmetric(cAdj)
            error('Graph must be symmetric');
        end
        
        zSel = zi > zj;

        zi = zi(zSel);
        zj = zj(zSel);
        zv = zv(zSel);
    end
    
    if opts.posOnly == 1
        fprintf('Dropping negative edges');
                
        zSel = zv > 0;

        zi = zi(zSel);
        zj = zj(zSel);
        zv = zv(zSel);
    end
    
    if opts.doFisher 
        fprintf('Fisher transform');
        zv = atanh(zv);
    end
    
    outNet = sprintf('%s/%s_net.txt',fDir,fName);
    outPart = sprintf('%s/%s_part.txt',fDir,fName);
    
    cTab.x = zi;
    cTab.y = zj;
    cTab.w = zv;
    cTab = struct2table(cTab);
    writetable(cTab,outNet,'WriteVariableNames',0,'delimiter','\t');
    
    outCmd = sprintf('java -jar %s mdl --verbose -g %s -o %s %s',opts.jarPath,outNet,outPart,opts.jarParam);
    
    [oStat,oRes] = system(outCmd);
    disp(oRes);
    
    inPart = fastTxtRead(outPart);
    
    cM.nodeIdx = str2double(inPart(:,1));
    cM.moduleIdx = str2double(inPart(:,2));
    assert(isequal(cM.nodeIdx,(1:G.numnodes)'));
    
    nodeNames = G.Nodes.Names;

    for i = 1:luniq(cM.moduleIdx)
        zSel = cM.moduleIdx+1 == i;
        outSetNames{i} = nodeNames(zSel);
    end
    
    if opts.doPrint
        disp(cellToColTable(outSetNames));
    end
    
    %%
    cAdjP = cAdj;
    cAdjP(cAdj < 0) = 0;
    cAdjN = abs(cAdj);
    cAdjN(cAdj > 0) = 0;
    
    mQ(2:3) = [ weightedModularityQ(cAdjP,cM.moduleIdx) -weightedModularityQ(cAdjN,cM.moduleIdx) ];
    mQ(1) = nansum(mQ);
    
%     if opts.clearTemp
%         rm(fDir)
%     end

end