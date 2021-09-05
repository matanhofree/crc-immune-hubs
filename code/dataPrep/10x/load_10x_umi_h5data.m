function [molData] = load_10x_umi_h5data(fileName,location)

    if nargin < 2
        location = '/';
    end
    dList = h5info(fileName,location);
    
    if isfield(dList,'Datasets') && ~isempty(dList.Datasets)
        colNames = { dList.Datasets.Name };    
        colNamesSt = matlab.lang.makeValidName(colNames);    
        for i = 1:length(colNames)

            molData.(colNamesSt{i}) = hdf5read(fileName,[ location colNames{i} ]);
        end
    end
    %%
    if isfield(dList,'Groups') && ~isempty(dList.Groups)                 
        
        groupNames = { dList.Groups.Name };
        groupNamesSt = matlab.lang.makeValidName(regexprep(groupNames,'/',''));
        for i = 1:length(groupNamesSt)
            molData.(groupNamesSt{i}) = load_10x_umi_h5data(fileName,[ groupNames{i} '/' ]);
        end
        
    end
    %%
    if isfield(dList,'Attributes') && ~isempty(dList.Attributes)
        attrNames = { dList.Attributes.Name };
        attrNamesSt = matlab.lang.makeValidName(attrNames);
        attrNamesSt = matlab.lang.makeUniqueStrings(attrNamesSt);
        
        attrNamesSt = regexprep(attrNamesSt,'(.*)','att_$1');
        
        for i = 1:length(attrNames)
            molData.(attrNamesSt{i}) = h5readatt(fileName,location,attrNames{i});
        end

    end

    
end


