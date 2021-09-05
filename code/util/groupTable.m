function outT = groupTable(cTable,groupVar,dropVar,inOpts)

    defaultOpts.maxCatVar = 10;
    defaultOpts.mergeNumeric = @(x)mean(x);
    defaultOpts.mergeCat = @(x)strjoin(x,';');
    defaultOpts.keyVar = 'sampleID';
    defaultOpts.outTable = 0;

    
    if (exist('inOpts','var') == 1)
        if isnumeric(inOpts) && inOpts == 1
            opts = defaultOpts;
            opts.outTable = 1;
        else
            opts = mergeOption(inOpts,defaultOpts);
        end
    else
        opts = defaultOpts;
    end
    clear defaultOpts % inOpts;


    if istable(cTable)
        cTable = table2struct(cTable,'ToScalar',1);
    end
% 
%     if ~exist('dropVar','var') || isempty(dropVar)
%         dropVar = 'sampleID';
%     end

    cField = fieldnames(cTable);
    
    if exist('dropVar','var') && ~isempty(dropVar)
        cField = setdiff(cField,dropVar);
    end
    
    if ischar(groupVar)
        if ~ismember(groupVar,cField)
            disp(cField)
            error('Group var is missing.')
        end
        groupV = cTable.(groupVar);
        cField = setdiff(cField,groupVar);
        opts.keyVar = groupVar;
    else
        if length(groupVar) ~= length(cTable.(cField{1}))
            error('Size missmatch');
        end
        groupV = groupVar;
    end
    
    
    % cTable = structSelectField(cTable,cField);

    [outT.(opts.keyVar),~,~,~,cPos] = fastUnique(groupV);

    mergeCat = opts.mergeCat;
    mergeNumeric = opts.mergeNumeric;

    for zi = 1:length(cField)
        cf = cField{zi};
        cVal = cTable.(cf);
        
        if iscategorical(cVal)
            cVal = cellstr(cVal);
        end
                    
        if iscell(cVal)
            outT.(cf) = cellfun(@(x)mergeCat(unique(cVal(x))),cPos,'unif',0);
        else
            outT.(cf) = cellfun(@(x)mergeNumeric(cVal(x)),cPos);
        end
    end

    if opts.outTable == 1
        outT = struct2table(outT);
    end
end