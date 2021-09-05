function [outS,outStubName] = load_10x_cr3_filteredBC_h5(dataFile,inOpts)


    defaultOpts.removeNA = 0;    
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end        
    clear defaultOpts;
    disp(opts);
    

    if ~isstruct(dataFile)        
        fprintf('Loading data from file: %s\n',dataFile);

        if exist(dataFile,'file')                     


            inData = load_10x_umi_h5data(dataFile);

            [~,outStubName] = fileparts(dataFile);
        else
            error('File not found');
        end
    else
        inData = dataFile;        
        outStubName = inputname(1);        
    end
    
    assert(inData.att_version == 2);
    %%
    outS.feature_type = categorical(get(inData.matrix.matrixfeatures.feature_type,'data'));
    outS.genome = categorical(get(inData.matrix.matrixfeatures.genome,'data'));
    outS.geneID = get(inData.matrix.matrixfeatures.name,'data');
    outS.ensgID = get(inData.matrix.matrixfeatures.id,'data');   
    
    
    outS.barcodes = get(inData.matrix.barcodes,'data');
    
    %%
    
    D = double(inData.matrix.shape(1));
    N = double(inData.matrix.shape(2));
    cNNZ = double(length(inData.matrix.data));    
    %%
    
    
    cDiff = diff(inData.matrix.indptr);
    %%
    cList = find(cDiff);
    cList = cList(:)';
    %%
    spColIdx = zeros(cNNZ,1);
    for ja = cList
       zidx = (inData.matrix.indptr(ja)+1):inData.matrix.indptr(ja+1);
       spColIdx(zidx) = ja;
       
       % assert(luniq(inData.matrix.indices(zidx)) == length(zidx),'Assumption about index uniquness is broken')
       % cRidx = inData.matrix.indices(jidx);
       % cRdata = inData.matrix.data(jidx);
       % adjustedCounts(cRidx,ja)
    end
    %% 
    
    outS.rawCount = sparse(double(inData.matrix.indices+1),spColIdx,double(inData.matrix.data),D,N);
    %%
    if isfield(inData,'att_chemistry_description')
        cv = inData.att_chemistry_description;
        outS.chemistry = repmat( { cv } ,N,1);
    end
    
    % cv = inData.att_library_ids;
%     if length(cv) == 1
%         cidx = str2double(regexprep(outS.barcodes,'[TCGAN]+-',''));
%         assert(unique(cidx) == 1);
%         outS.barcodes = regexprep(outS.barcodes,'([^-]+)-.*','$1');
%         outS.sampleID = mergeStringPair('%s_id-%s',cv,outS.barcodes);
%     else
%         error('Not yet implemented');
%     end
%     outS.ambient_expression = inData.matrix.ambient_expression;
%     outS.latent_scale = inData.matrix.latent_scale;
%     outS.latent_cell_probability = inData.matrix.latent_cell_probability;
%     %%
    
    
end