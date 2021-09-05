function optionFinal = mergeOption(option,optionDefault,mergeDefaultEmpty)
% Merge two struct options into one struct
% Usage:
% optionFinal = mergeOption(option,optionDefault)
    if nargin < 3
        mergeDefaultEmpty = 0;
    end
    
    if mergeDefaultEmpty == 0
        optionFinal=optionDefault;
    end
    if isempty(option)
        return;
    end

    names=fieldnames(option);
    for i=1:numel(names)
        cfield = option.(names{i});
        if isstring(cfield)
            cfield = cellstr(cfield);
        end
        
        if isfield(optionDefault,names{i})
            if isstruct(cfield)
                cfield = mergeOption(cfield,optionDefault.(names{i}),mergeDefaultEmpty);            
%                 if isstring(cfield)
%                     cfield = cellstr(cfield);
%                 end                

                optionFinal.(names{i}) = cfield;
            elseif mergeDefaultEmpty == 1 && isempty(cfield)                
                optionFinal.(names{i}) = optionDefault.(names{i});
            else
                optionFinal.(names{i}) = cfield;
            end
        else            
            optionFinal.(names{i}) = cfield;
        end
    end

end