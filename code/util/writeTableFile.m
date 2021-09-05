function cTab = writeTableFile(inMat,outFile,inRow,inCol,inOpts)

    defaultOpts.fixRow = 0;
    defaultOpts.fixCol = 0;
    defaultOpts.uniqueRow = 0;
    defaultOpts.uniqueCol = 1; 
    defaultOpts.forceRow = 0; 
    defaultOpts.forceCol = 0; 
    defaultOpts.writeRowNames = 1;
    defaultOpts.writeColNames = 1;
    defaultOpts.makeDir = 1;    
    defaultOpts.forceClobber = 0;
    defaultOpts.suffixDelim = 1;
    defaultOpts.delim = '\t';
    defaultOpts.quoteStrings = false;
    defaultOpts.fileType = 'text';

        
    if (exist('inOpts','var') == 1)
        if isstruct(inOpts)
            opts = mergeOption(inOpts,defaultOpts);
        elseif ismatrix(inOpts)
            opts = defaultOpts;
            opts.forceClobber = inOpts;
        end            
    else
        opts = defaultOpts;
    end
    clear defaultOpts
    
    if ~exist('outFile','var') || isempty(outFile)
        outFile = sprintf('./temp/%s.tsv',inputname(1));
        opts.forceClobber = 1;
    end
    
    if istable(inMat) && ~isempty(inMat.Properties.RowNames)
        inRow = inMat.Properties.RowNames;
    end

    writeCol = opts.writeColNames;
    writeRow = opts.writeRowNames;                   
    
    if ~exist('inRow','var') || isempty(inRow)
        inRow = [];        
        if opts.forceRow == 0
            writeRow = 0;
        end            
    elseif isnumeric(inRow) && length(inRow) == 1
        opts.forceClobber = 1;
        inRow = [];
        
    end        
    
    if ~exist('inCol','var') || isempty(inCol)
        inCol = [];
        
        if istable(inMat)
            inCol = inMat.Properties.VariableNames;
        else             
            if opts.forceCol == 0
                writeCol = 0;
            end                
        end
    end
    
    colID = inCol;
    rowID = inRow;
    
    if ~exist('outFile','var')
        outFile = [];        
    end
    
    outFileRename = [];
    if ~isempty(outFile)        
        [baseDir,fname,ext] = fileparts(outFile);
        
        if ~exist(baseDir,'dir') && opts.makeDir             
            p = mkdir(baseDir);            
            if p ~= 1
                error('Unable to create dir');
            end            
        end                
        
        if isempty(ext) 
            if exist([ baseDir filesep fname],'dir') 
                
                baseDir = [ baseDir filesep fname];
                fname = inputname(1);
                ext = '.csv'
                if ~isempty(fname)
                    error('Input file name required');
                end
            else
                ext = '.csv';
            end
        elseif ismember(ext,{'.xls', '.xlsx'})
           ext = '.xlsx';
           opts.fileType = 'spreadsheet';
            
        elseif ~ismember(ext,{ '.csv','.txt','.dat' })
            outFileRename = sprintf('%s/%s%s',baseDir,fname,ext);
            ext = '.txt';
        
        end
        outFile = sprintf('%s/%s%s',baseDir,fname,ext);
    end
    
    [N,D] = size(inMat);                                  

    if writeRow        
        if isempty(rowID)            
            rowID = mergeStringPair('r',1:N);
        else
            if opts.uniqueRow && ~isunique(rowID)
                rowID = matlab.lang.makeUniqueStrings(rowID);
            end
        end

    end
                
    if writeCol         
        if isempty(colID)            
            colID = mergeStringPair('c',1:D);
        else
            if opts.uniqueCol && ~isunique(colID)
                colID = matlab.lang.makeUniqueStrings(colID);
            end
        end
    end
    
     if writeCol && length(colID) ~= D 
        error('Column length missmatch');
    end

    if writeRow && length(rowID) ~= N 
        error('Row length missmatch');
    end
    
    if  writeCol && opts.fixCol && ~isempty(colID)
        colID = matlab.lang.makeValidName(colID);
        if ~isequal(colID,inCol)
            fprintf('Fixed column names\n')
        end
    end
    
    if writeRow && opts.fixRow  && ~isempty(rowID)
        rowID =  matlab.lang.makeValidName(rowID);
        if ~isequal(rowID,inRow)
            fprintf('Fixed row names\n');
        end        
    end

    %% Create Table
    if isstruct(inMat)

        if ~isempty(colID)
            cTab = struct2table(inMat,'VariableNames',colID);                                    
        else
            cTab = struct2table(inMat);
            writeCol = opts.writeColNames;
        end   
    elseif istable(inMat)       
            
        % if ~isempty(colID)
        %     cTab.Proporties.VariableNames = colID;                                                                 
        % end
        % writeRow = 1;
        cTab = inMat;
    elseif ismatrix(inMat)
        if iscell(inMat)        
            if ~isempty(colID)
                cTab = cell2table(inMat,'VariableNames',colID);                                    
            else
                cTab = cell2table(inMat);                                    
            end
        else
            if ~isempty(colID)
                cTab = array2table(inMat,'VariableNames',colID);                                    
            else
                cTab = array2table(inMat);                                    
            end        
        end

    
    else
        error('Unable to write inMat')
    end    
    
    if ~istable(inMat) && writeRow && ~isempty(rowID) 
        if isunique(rowID) 
            cTab.Properties.RowNames = rowID(:);
        else
            cTab = [ cell2table(rowID,'VariableNames',{'id'}) cTab ];
            writeRow = 0;
        end
    end
    
    if ~exist(outFile,'file') || opts.forceClobber == 1
        % cTab(1:10,1:10)
        if strcmp(opts.fileType,'spreadsheet')
            writetable(cTab,outFile,'WriteRowNames',writeRow,'WriteVariableNames',writeCol,'FileType',opts.fileType);
        else
            writetable(cTab,outFile,'WriteRowNames',writeRow,'WriteVariableNames',writeCol,'delimiter',opts.delim,'QuoteStrings',opts.quoteStrings,'FileType',opts.fileType);
        end
        if ~isempty(outFileRename)
            movefile(outFile,outFileRename);
        end
    else
        warning('File %s already exits. Remove file, change write location, or set forceClobber to true to write\n',outFile);
    end        
    
    
    if nargout == 0
        cTab = [];
    end
    
end