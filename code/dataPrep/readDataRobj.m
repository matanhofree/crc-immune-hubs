function outData = readDataRobj(dirPath)


    defaultOpts.nRow = 'geneID';
    defaultOpts.nCol = 'sampleID';

    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;


    [fList,~,fName] = fileList(dirPath);

    varTypeList = { '_dZ_' '_dSp_' '_dDvec' '_dDmat' };

    cSel = strgrep(fName,'_dSp_');
    if sum(cSel) > 0
        fList = [ fList(~cSel); fList(cSel) ];
        fName = [ fName(~cSel); fName(cSel) ];
    end

    for zi = 1:length(fList)

        cf = fList{zi};
        cname = fName{zi};

        if ~strgrepE(cname,varTypeList)
            fprintf('Warn: unknown variabl type. Skipping %s\n',cf);
        else
            fprintf('Reading: %s\n',cname);
        end


        if strgrep(cname,'_dSp_')

            cStub = regexprep(cname,'.*_dSp_([^.]*)\..*','$1');
            cMat = load_h5data(cf);

            if all(ismember({'i','i','v'},fieldnames(cMat)))
                fprintf('Processing %s as matlab h5 sparse matrix\n',cStub);
                if all(ismember({'nRow','nCol'},fieldnames(cMat)))
                    outData.(cStub) = sparse(cMat.i,cMat.j,cMat.v,cMat.nRow,cMat.nCol);
                elseif ~isempty(opts.nRow) && ~isempty(opts.nCol)
                    if ischar(opts.nRow)
                        opts.nRow = length(outData.(opts.nRow));
                    end
                    if ischar(opts.nCol)
                        opts.nCol = length(outData.(opts.nCol));
                    end

                    outData.(cStub) = sparse(cMat.i,cMat.j,cMat.v,opts.nRow,opts.nCol);
                else
                    fprintf('Warn: sparse variable does not specify dimensions\n');
                    outData.(cStub) = sparse(cMat.i,cMat.j,cMat.v);
                end
            else
                fprintf('Attempting to process %s as cellranger file\n',cStub);
                error('Not yet implemented');
                cMat = load_10x_cr3_filteredBC_h5_T(cMat);
            end

            fprintf('Read %s (%dx%d) sparse matrix\n',cStub,size(outData.(cStub),1),size(outData.(cStub),2));
        end

        if strgrep(cname,'_dZ_')

            cStub = regexprep(cname,'.*_dZ_([^.]*)\..*','$1');
            isZip = strgrep(cf,'.gz$');

            if isZip
                gunzip(cf)
                cf = regexprep(cf,'.gz$','');
            end
            cMat = readtable(cf);
            if isZip
                gzip(cf);
            end
        
            fprintf('Read %s (%dx%d) table\n',cStub,size(cMat,1),size(cMat,2));
    
            outData.(cStub) = table2struct(cMat,'toScalar',1);
        end


        if strgrep(cname,'_dDvec_')

            cStub = regexprep(cname,'.*_dDvec_([^.]*)\..*','$1');
            
            cMat = fastTxtRead(cf);
            
            fprintf('Read %s (%dx%d) variable\n',cStub,size(cMat,1),size(cMat,2));
    
            outData.(cStub) = cMat;
        end

        if strgrep(cname,'_dDmat_')

            cStub = regexprep(cname,'.*_dDmat_([^.]*)\..*','$1');
            
            [cMat,cErr,cRowH,colH] = fastMatRead(cf,',');
            
            fprintf('Read %s (%dx%d) numeric matrix\n',cStub,size(cMat,1),size(cMat,2));
    
            outData.(cStub) = cMat;
        end


    end

end