function leafOrder = constrainedLeafOrder(expSub,clustV,typeOrd,cDist,cLink,useLeafOrd,inOpts)

    defaultOpts.distPower = 1;
    defaultOpts.sortClustNames = 1;

    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;

    if isempty(clustV)
        fprintf('Unconstrained clustering.');

        cPDist = pdist(expSub,cDist);
        cLinkM = linkage(cPDist,cLink);

        if useLeafOrd
            leafOrder = optimalleaforder(cLinkM,cPDist);
        else
            % [~,~,leafOrder] = dendrogram(cLinkM);
            [~,cOrder] = dendrogram_order(cLinkM);
        end

    else
        cN = size(expSub,1);
        assert(cN == length(clustV));

        [clist,~,~,~,clPos] = fastUnique(clustV);
        if opts.sortClustNames
            [clist,cidx] = sort(clist);
            clPos = clPos(cidx);
        end

        if typeOrd == 2
            fprintf('Block-wise HC \n');

            blockExp = cell2mat(cellfun(@(x)mean(expSub(x,:),1),clPos,'unif',0));
            cPDist = pdist(blockExp,cDist);


            if opts.distPower > 1
                cPDist = cPDist.^opts.distPower;
            end

            cLinkM = linkage(cPDist,cLink);
            cOrder = optimalleaforder(cLinkM,cPDist);
%             if useLeafOrd
%                 cOrder = optimalleaforder(cLinkM,cPDist);
%             else
%                 % [~,~,cOrder] = dendrogram(cLinkM);
%                 [~,cOrder] = dendrogram_order(cLinkM);
% 
%             end
            clist = clist(cOrder);
            clPos = clPos(cOrder);
        end


        leafOrder = zeros(cN,1);
        clN = length(clist);

        fprintf('Per block optimal leaf order\n');
        zp = 1;
        for i = 1:clN
            fprintf('Running %d of %d (%d smp)\n',i,clN,length(clPos{i}));
            ci = clPos{i};

            nC = length(ci);

            if nC > 1

                cPDist = pdist(expSub(ci,:),cDist);


                if opts.distPower > 1
                    cPDist = cPDist.^opts.distPower;
                end

                cLinkM = linkage(cPDist,cLink);

                if useLeafOrd
                    cOrder = optimalleaforder(cLinkM,cPDist);
                else
                    % [~,~,cOrder] = dendrogram(cLinkM);
                    [~,cOrder] = dendrogram_order(cLinkM);

                end
                leafOrder(zp:zp+nC-1) = ci(cOrder);
            else
                leafOrder(zp:zp+nC-1) = ci;
            end

            zp = zp+nC;
        end


    end
end

%%
%         if typeOrd == 1
%             fprintf('Per block optimal leaf order\n');
%
%             zp = 1;
%             for i = 1:clN
%                 fprintf('Running %d of %d (%d smp)\n',i,clN,length(clPos{i}));
%                 ci = clPos{i};
%
%                 nC = length(ci);
%
%                 if nC > 1
%
%                     cPDist = pdist(expSub(ci,:),cDist);
%
%
%                     if opts.distPower > 1
%                         cPDist = cPDist.^opts.distPower;
%                     end
%
%                     cLinkM = linkage(cPDist,cLink);
%
%                     if useLeafOrd
%                         cOrder = optimalleaforder(cLinkM,cPDist);
%                     else
%                         % [~,~,cOrder] = dendrogram(cLinkM);
%                         [~,cOrder] = dendrogram_order(cLinkM);
%
%                     end
%                     leafOrder(zp:zp+nC-1) = ci(cOrder);
%                 else
%                     leafOrder(zp:zp+nC-1) = ci;
%                 end
%
%                 zp = zp+nC;
%             end
%         elseif typeOrd == 2
%             fprintf('Per block optimal leaf order with HC on top PCs\n');
%
%             zp = 1;
%             for i = 1:clN
%                 fprintf('Running %d of %d (%d smp)\n',i,clN,length(clPos{i}));
%                 ci = clPos{i};
%
%                 nC = length(ci);
%
%                 if nC > 1
%
%                     cPDist = pdist(expSub(ci,:),cDist);
%
%
%                     if opts.distPower > 1
%                         cPDist = cPDist.^opts.distPower;
%                     end
%
%                     cLinkM = linkage(cPDist,cLink);
%
%                     if useLeafOrd
%                         cOrder = optimalleaforder(cLinkM,cPDist);
%                     else
%                         % [~,~,cOrder] = dendrogram(cLinkM);
%                         [~,cOrder] = dendrogram_order(cLinkM);
%
%                     end
%                     leafOrder(zp:zp+nC-1) = ci(cOrder);
%                 else
%                     leafOrder(zp:zp+nC-1) = ci;
%                 end
%
%                 zp = zp+nC;
%             end
%         else
%             error('Not implemented');
% %             leafOrderX = 1:size(expSub,2);
%         end
%%%

%             cPDist = squareform(pdist(expSub',cDist));
%             zfact = 10;
%             zPDistReNorm = cPDist*zfact;
%
%             %%
%             [cCLpos.uniqueList,~,~,~,cCLpos.pos] = fastUnique(clustV);
%             for i = 1:length(cCLpos.uniqueList)
%                 zc = cCLpos.pos{i};
%                 zPDistReNorm(zc,zc) = cPDist(zc,zc);
%             end
%             %%
%
%             zPDistReNorm = squareform(zPDistReNorm);
%
%             cLinkM = linkage(zPDistReNorm,cLink);
%             %%
%
%             leafOrderX = optimalleaforder(cLinkM,zPDistReNorm);
%
%%
%        elseif typeOrd == 2
%             fprintf('Using k cluster opimal order heuristic\n');
%
%
%             [cMean,~,cBatch] = summarize_subset_value(expSub,clustV);
%
%             %%
%             cPDist = (pdist(cMean',cDist));
%             cLinkM = linkage(cPDist,'complete');
%
%             leafOrderX = optimalleaforder(cLinkM,cPDist);
%
%             lBatchOrder = cBatch(leafOrderX)
%             %%
%             [cCLpos.uniqueList,~,~,~,cCLpos.pos] = fastUnique(clustV);
%             [~,zia,zib] = intersect(lBatchOrder,cCLpos.uniqueList);
%             %
%             reOrd(zia) = zib;
%             %
%             assert(isequal(cCLpos.uniqueList(reOrd),lBatchOrder));
%             %
%             cCLpos.uniqueList = cCLpos.uniqueList(reOrd);
%             cCLpos.pos = cCLpos.pos(reOrd);
%             %
%
%             for i = 1:length(cCLpos.uniqueList)
%                 zc = cCLpos.pos{i};
%                 cPDist = pdist(expSub(:,zc)',cDist);
%
%                 cLinkM = linkage(cPDist,cLink);
%                 leafOrderX = optimalleaforder(cLinkM,cPDist);
%
%                 globalReorder{i} = zc(leafOrderX);
%             end
%             %%
%             leafOrderX = [ globalReorder{:} ];