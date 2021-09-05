function strOut = structSubSelectMat(strIn,selectV,applyToCell,depthR)
% Apply a selection on all variables in a struct which match a particular     

    strOut = strIn;   
    
    if isempty(selectV) 
        fprintf('Selection is empty\n');
        return;
    end
    
    if iscellstr(selectV) && isfield(strIn,'sampleID')
        selSub = ismember(strIn.sampleID,selectV);

        if any(selSub)
            fprintf('Using sampleID as selection dim (%d/%d)\n',sum(selSub),length(selSub));
            selectV = selSub;
        end
    end
    
    if iscellstr(selectV) && isfield(strIn,'geneID')
        selSub = ismember(strIn.geneID,selectV);

        if any(selSub)
            fprintf('Using geneID as selection dim (%d/%d)\n',sum(selSub),length(selSub));
            selectV = selSub;
        end
    end
    
    if iscellstr(selectV)
        error('Improper selection criteria\n');
    end
    
    if all(selectV)
        fprintf('Nothing to subset (selecting all)\n')
        return
    end
    

    if sum(selectV) == 0
        error('Selection may be removing all');
    end
    
    if nargin < 3
%         if length(selectV) < 100 
%             applyToCell = 0;
%         else
%             applyToCell = 1;
%         end
        applyToCell = 0;
    end
    
    if nargin < 4
        depthR = inf;
    end
    
    if depthR < 1
        strOut = strIn;
        return;
    end
    
    lenV = length(selectV);

    fnames = [];
    ctype = 0;
%     if applyToCell
%         fprintf('Filtering cells\n');
%     end   
    if iscell(strIn) 
        ctype = 1;
        nEl = numel(strIn);    
    elseif isstruct(strIn) 
        ctype = 2;
        fnames = fieldnames(strIn);
        nEl = length(fnames);
    elseif ismatrix(strIn)
        ctype = 3;
        nEl = 1;
    else
        error('Type unknown');
    end
    
    for i = 1:nEl
        
        if ctype == 2 
            cname = fnames{i};
            zmat = strIn.(cname);                        
        elseif ctype == 1
            zmat = strIn{i};
        elseif ctype == 3
            zmat = strIn;
        end

        if (isstruct(zmat))
            zmat = structSubSelectMat(zmat,selectV,applyToCell,depthR-1);
        else
            switch ndims(zmat)
                case 1
                    if (length(zmat) == lenV)
                        % zmat(~selectV) = [];
                        zmat = zmat(selectV);
                    elseif iscell(zmat) && applyToCell                        
                        zmat = cellfun(@(subMat)structSubSelectMat(subMat,selectV,applyToCell,depthR-1),zmat,'uniformoutput',0);
                    end
                case 2
                    if (size(zmat,1) == lenV)
                        % zmat(~selectV,:) = [];
                        zmat = zmat(selectV,:);
                    end
                    if (size(zmat,2) == lenV)
                        % zmat(:,~selectV) = [];
                        zmat = zmat(:,selectV);
                    end

                    if iscell(zmat) && applyToCell
                        zmat = cellfun(@(subMat)structSubSelectMat(subMat,selectV,applyToCell,depthR-1),zmat,'uniformoutput',0);
                    end                    
                case 3
                    if (size(zmat,1) == lenV)
                        % zmat(~selectV,:,:) = [];
                        zmat = zmat(selectV,:,:);
                    end
                    if (size(zmat,2) == lenV)
                        % zmat(:,~selectV,:) = [];
                        zmat = zmat(:,selectV,:);
                    end
                    if (size(zmat,3) == lenV)
                        % zmat(:,:,~selectV) = [];
                        zmat = zmat(:,:,selectV);
                    end

                    if iscell(zmat) && applyToCell                        
                        zmat = cellfun(@(subMat)structSubSelectMat(subMat,selectV,applyToCell,depthR-1),zmat,'uniformoutput',0);
                    end

                otherwise
                    if iscell(zmat) && applyToCell                        
                        zmat = cellfun(@(subMat)structSubSelectMat(subMat,selectV,applyToCell,depthR-1),zmat,'uniformoutput',0);
                    else
                        warning('Might be error - skipping something - double check');
                        continue;
                    end 
            end
        end
        
        if ctype == 2 
            strOut.(cname) = zmat;        
        elseif ctype == 1
            strOut{i} = zmat;        
        elseif ctype == 3
            strOut = zmat;
        end
        
    end
end