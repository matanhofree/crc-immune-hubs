function wManifest = writeDataToR(inData,outDir,inName,skipVars,inOpts)
    
    defaultOpts.defName = 'data';
    defaultOpts.gzip = 1;

    defaultOpts.structSep = '_dZ_';
    defaultOpts.sparseSep = '_dSp_';
    defaultOpts.denseVectorSep = '_dDvec_';
    defaultOpts.denseMatSep = '_dDmat_';
    
    defaultOpts.charSep = '_chr_';   
    defaultOpts.writeSparse = 'HDF5';
    defaultOpts.writeDense = 'HDF5';

    defaultOpts.removeZeros = 1;
    defaultOpts.inSub = 0;
    
    defaultOpts.skipData = 0;

    defaultOpts.delim = ',';
    defaultOpts.suffix = '.csv';
    
    defaultOpts.clobber = 0;
    
    
                
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    
    if ~exist('skipVars','var') || isempty(skipVars)
        skipVars = {};
    elseif ischar(skipVars)
        skipVars = strsplit(skipVars,',');        
    end
    
    if ~exist('inName','var') || isempty(inName)
        inName = inputname(1);     
    end
    
    if isempty(inName)       
        inName = opts.defName;        
    end           
            
    if ~exist(outDir,'dir') 
        [~,~,~] = mkdir(outDir);
    end
    
    if isempty(opts.suffix)
        if strcmp(opts.delim,'\t')
            opts.suffix = '.tsv';
        else 
            opts.suffix = '.csv';
        end
    end
    wManifest = [];
    fnames = fieldnames(inData);      

%     wManifest.varName
%     wManifest.varDim
%     wManifest.varType         
    zx = 1;
    for i = 1:length(fnames)
        
        cname = fnames{i};   
        
        if ismember(cname,skipVars)            
            fprintf('Skipping: %s (inList)\n',cname);
            continue;
        end
                
        cMat = inData.(cname);
        [zD,zN] = size(cMat);                                        
        
        if isstruct(cMat)                         
            
            cLen = structfun(@(x)size(x,1),cMat);
            subName = [ outDir '/' inName opts.structSep cname ];
            isTable = 0;
            
            if all(cLen == cLen(1))                                
                try 
                    cTable = struct2table(cMat);  
                    
                    if opts.clobber == 0 && (exist([ subName opts.suffix ],'file') || exist([ subName opts.suffix '.gz' ],'file'))
                        fprintf('Skipping %s. File exists %s\n',cname,[ subName opts.suffix ]);
                        continue;
                    end                    
                    
                    writetable(cTable,[ subName opts.suffix ],'Delimiter',opts.delim);                    
                                        
                    if opts.gzip
                        gzip([ subName opts.suffix ]);
                        delete([ subName opts.suffix ]);
                    end
                    
                    isTable = 1;
                    
                    wManifest.varName{zx,1} = subName;
                    wManifest.varDim{zx,1} = size(cTable);
                    wManifest.varType{zx,1} = 'table';
                    zx = zx + 1;
                catch e
                    fprintf('Struct %s is not a simple table. Attempting to write as plain struct.\n',cname);
                end                
                
                if isTable
                    continue;
                end
            end
            
            opts.inSub = 1;
            subManif = writeDataToR(cMat,subName,cname,'sampleID',opts);                        
            
            wManifest.varName{zx,1} = subName;
            wManifest.varDim{zx,1} = [ NaN NaN ];
            wManifest.varType{zx,1} = 'struct';
            zx = zx + 1; 
            
            if ~isempty(subManif) && length(subManif.varName) > 0
                wManifest.varName = [wManifest.varName; subManif.varName];
                wManifest.varDim = [wManifest.varDim; subManif.varDim];
                wManifest.varType = [wManifest.varType; subManif.varType];

                zx = zx + length(subManif.varName);
            end
            
        elseif issparse(cMat) 
            
            subName = [ outDir '/' inName opts.sparseSep cname ];
            
            if opts.skipData 
                fprintf('Skipping writing of data matrix (%s)\n',cname);                    
                continue;
            end
            
                                
            if opts.clobber == 0 && exist([ subName '.h5' ],'file')
                fprintf('Skipping %s. File exists %s\n',cname,[ subName '.h5' ]);
                continue;
            end                    
     
            
            switch opts.writeSparse                
                case 'HDF5'
                    fprintf('Writing %s in hdf5 (matlab) format\n',cname);                    
                    outS = [];
                    [ outS.i, outS.j, outS.v ] = find(cMat);
                    % outS.nVar = nnz(cMat);
                    [outS.nRow,outS.nCol] = size(cMat);
                    save([ subName '.h5'],'-v7.3','-struct','outS');                    
                    
                    wManifest.varName{zx,1} = subName;
                    wManifest.varDim{zx,1} = [ zD zN ]
                    wManifest.varType{zx,1} = 'sparse-h5';   
                    zx = zx + 1; 

                otherwise
                    fprintf('Skipping: %s. Choose a sparse format to use\n',cname);
                    
            end
    
        else
            % Not sparse            
  
            if min([zD zN]) == 1 % Vector               
                subName = [ outDir '/' inName opts.denseVectorSep cname ];

                if opts.clobber == 0 && (exist([ subName opts.suffix ],'file') || exist([ subName opts.suffix '.gz' ],'file'))
                    fprintf('Skipping %s. File exists %s\n',cname,[ subName opts.suffix ]);
                    continue;
                end   
                
                try 
                    cTable = array2table(cMat);                    
                    writetable(cTable,[ subName opts.suffix ],'Delimiter',opts.delim,'WriteVariableNames',0);                    
                                        
                    if opts.gzip
                        gzip([ subName opts.suffix ]);
                        delete([ subName opts.suffix ]);
                    end
                                        
                    wManifest.varName{zx,1} = subName;
                    wManifest.varDim{zx,1} = size(cTable);
                    if iscell(cMat)
                        wManifest.varType{zx,1} = 'vector-char';
                    elseif isnumeric(cMat)
                        wManifest.varType{zx,1} = 'vector-num';
                    else
                        wManifest.varType{zx,1} = 'vector-other';
                    end                                        
                    zx = zx + 1; 
    
                catch e
                    fprintf('Vector %s could not be written - skipping.\n',cname);
                end 
                
            else % Matrix
                
                subName = [ outDir '/' inName opts.denseMatSep cname ];

                if opts.clobber == 0 && (exist([ subName opts.suffix ],'file') || exist([ subName opts.suffix '.gz' ],'file') || exist([ subName '.h5' ],'file'))
                    fprintf('Skipping %s. File exists %s\n',cname,[ subName ]);
                    continue;
                end   
                
                wManifest.varName{zx,1} = subName;                                        
                wManifest.varDim{zx,1} = size(cMat);
                
                if iscell(cMat)
                    
                    cTable = cell2table(cMat);
                    
                    writetable(cTable,[ subName opts.suffix ],'Delimiter',opts.delim,'WriteVariableNames',0);                    

                    if opts.gzip
                        gzip([ subName opts.suffix ]);
                        delete([ subName opts.suffix ]);
                    end
                    wManifest.varType{zx,1} = 'matrix-char';                                    
                    
                elseif isnumeric(cMat)
                    
                    wManifest.varType{zx,1} = 'matrix-num-h5';
                    
                    outS = [];
                    outS.data = cMat;
                    save([ subName '.h5'],'-v7.3','-struct','outS');                    
                    
                else                  
                    wManifest.varType{zx,1} = 'matrix-other';
                    
                    outS = [];
                    outS.data = cMat;
                    save([ subName '.h5'],'-v7.3','-struct','outS');                    

                end  
                zx = zx + 1; 

            end
        end
    end 

end


