function strOut = structSortMat(strIn,selectV,applyToCell)
% Apply a selection on all variables in a struct which match a particular     

    if nargin < 3
        applyToCell = 0;
    end

    fnames = fieldnames(strIn);
    lenV = length(selectV);    
    strOut = strIn;
    
    dropnan = @(x)x(~isnan(x));
    
    for i = 1:length(fnames)
        cname = fnames{i};
        
        zmat = strIn.(cname);                
        
        if (isstruct(zmat))
            zmat = structSortMat(zmat,selectV);
        else
            switch ndims(zmat)
                case 1
                    if (length(zmat) == lenV)
                        zmat = zmat(dropnan(selectV));
                    end
                case 2
                    if (size(zmat,1) == lenV)
                        zmat = zmat(dropnan(selectV),:);                    
                    end
                    if (size(zmat,2) == lenV)
                        zmat = zmat(:,dropnan(selectV));
                    end
                case 3
                    if (size(zmat,1) == lenV)
                        zmat = zmat(dropnan(selectV),:);
                    end
                    if (size(zmat,2) == lenV)
                        zmat = zmat(:,dropnan(selectV),:);
                    end
                    if (size(zmat,3) == lenV)
                        zmat = zmat(:,:,dropnan(selectV));
                    end                    
                otherwise
                    if iscell(zmat) && applyToCell                        
                        zmat = cellfun(@(subMat)structSortMat(subMat,selectV,applyToCell),zmat);
                    else
                        continue;
                    end                    
            end
        end
        strOut.(cname) = zmat;
    end

end