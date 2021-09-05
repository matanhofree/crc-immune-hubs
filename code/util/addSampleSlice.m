function [outS,outSlice] = addSampleSlice(inS,subSlice,sliceVar,sliceInS,mergeIn)

    if ~isstruct(inS) || ~isstruct(subSlice)
       error('Both inputs must be structs');
    end
    
    if ~exist('sliceVar','var') || isempty(sliceVar)
        sliceVar = 'sampleID';
    end
    
    if ~exist('sliceInS','var') || isempty(sliceInS)
        sliceInS = 0;
    end
    
    if ~exist('mergeIn','var') || isempty(mergeIn)
        mergeIn = 0;
    end
    
    if ~isfield(inS,sliceVar) || ~isfield(subSlice,sliceVar)
       error('Both inputs must be structs with the field - %s included',sliceVar);
    end
       
    idxSub = ismember(inS.(sliceVar),subSlice.(sliceVar));
    if sum(idxSub) ~= length(inS.(sliceVar))
        if sum(idxSub) == 0
            error('No overlap on sliceVar.\n');
        end
        
        if sliceInS == 0
            error('The slice struct (variable 2), is missing some ids from inS (variable 1).\n');
            
        elseif sliceInS == -1                       
            fprintf('Sub slicing variable 1 to match 2\n');    
            outS = structSubSelectMat(idxSub);
            
        elseif sliceInS == 1                       
            
            fprintf('Attempting to expand to match2\n');    
            
            [~,zia,zib] = intersect(inS.(sliceVar),subSlice.(sliceVar));
            
            cField = fieldnames(subSlice);
            cField(ismember(cField,sliceVar)) = [];
            subExp.(sliceVar) = inS.(sliceVar);
            
            N = length(inS.(sliceVar));
            
            for i = 1:length(cField)
                
                cf = cField{i};            
                cVal = subSlice.(cf);
                
                if iscell(cVal)
                    % zf = cell(N,1);
                    zf = repmat({ 'NA' },N,1);                    
                    zf(zia) = cVal(zib);
                else
                    zf = nan(N,1);       
                    zf(zia) = cVal(zib);                    
                end
                
                subExp.(cf) = zf;
            end        
            
            subSlice = subExp;
            outS = inS; 
        end            
    else                
        outS = inS;    
    end
    
        
    idxSub = ismember(subSlice.(sliceVar),outS.(sliceVar));    
    smpA = subSlice.(sliceVar)(idxSub);   
    %%
    if ~isempty(setdiff(outS.(sliceVar),smpA))
        error('Varible ids missmatch. Perhaps set sliceInS to 1\n');
    end
    
    subSlice = structSubSelectMat(subSlice,idxSub);    
    if ~isequal(outS.(sliceVar),subSlice.(sliceVar))
        [~,ia,ib] = intersect(outS.(sliceVar),subSlice.(sliceVar));
        subSlice = structSortMat(subSlice,ib(argsort(ia)));
    end
    assert(isequal(outS.(sliceVar)(:),subSlice.(sliceVar)(:)),'Funny buisness these ids should now match')
    
    subSlice = rmfield(subSlice,sliceVar);
    %%
    subSliceName = inputname(2);
    if isempty(subSliceName)
        subSliceName = 'subSlice';
    end
    
    if mergeIn == 1 || mergeIn == -1
        if mergeIn == 1
            bothF = intersect(fieldnames(inS),fieldnames(subSlice));
            if ~isempty(bothF)
                subSlice = rmfield(subSlice,bothF{:});        
            end
        end    
        
        flist = fieldnames(subSlice);
        for i = 1:length(flist)
            outS.(flist{i}) = subSlice.(flist{i});
        end        
    else
        outS.(subSliceName) = subSlice;
    end
    
end

