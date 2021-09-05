classdef discColorMapper
    properties
        nameSpace
        nameSpaceMap
        colorSpace        
        opts
        
    end
    
    methods
        function obj = discColorMapper(strList,colorSpaceIn,inOpts)
        
            defaultOpts.baseColor = [ 0.63 0.63 0.63 ];
            defaultOpts.sortNames = [];
            defaultOpts.errorOnNew = 0;
            defaultOpts.wrapColors = 1;

            if (exist('inOpts','var') == 1)
                opts = mergeOption(inOpts,defaultOpts);
            else
                opts = defaultOpts;
            end           
            
            % Prepare name space 
            if isempty(opts.sortNames)
                if isunique(strList)
                    nameSpace = strList;
                else
                    nameSpace = fastUnique(strList);
                end
            else 
                if opts.sortNames
                    nameSpace = unique(strList)
                else
                    nameSpace = fastUnique(strList)
                end                    
            end
            
            nSpace = length(nameSpace);
            
            if ~exist('colorSpaceIn','var') || isempty(colorSpaceIn)
                zColSet.rainbow14={'#882E72', '#B178A6', '#D6C1DE', '#1965B0', '#5289C7', '#7BAFDE', '#4EB265', '#90C987', '#CAE0AB', '#F7EE55', '#F6C141', '#F1932D', '#E8601C', '#DC050C'};
                zColSet.rainbow15={'#114477', '#4477AA', '#77AADD', '#117755', '#44AA88', '#99CCBB', '#777711', '#AAAA44', '#DDDD77', '#771111', '#AA4444', '#DD7777', '#771144', '#AA4477', '#DD77AA'};
                zColSet.rainbow18={'#771155', '#AA4488', '#CC99BB', '#114477', '#4477AA', '#77AADD', '#117777', '#44AAAA', '#77CCCC', '#777711', '#AAAA44', '#DDDD77', '#774411', '#AA7744', '#DDAA77', '#771122', '#AA4455', '#DD7788'};
                zColSet.rainbow21={'#771155', '#AA4488', '#CC99BB', '#114477', '#4477AA', '#77AADD', '#117777', '#44AAAA', '#77CCCC', '#117744', '#44AA77', '#88CCAA', '#777711', '#AAAA44', '#DDDD77', '#774411', '#AA7744', '#DDAA77', '#771122', '#AA4455', '#DD7788'};
                zColSet.crcTypeColor = { '#E41A1C' '#FB9A99' '#0072B2' '#B2DF8A' };
                colorSet = @(zCol)cell2mat(cellfun(@(x)hexColor(x(2:end)),zCol,'uniformoutput',0))';

                if nSpace <= 4
                    colorSpaceIn = colorSet(zColSet.crcTypeColor);
                elseif nSpace >= 5 && nSpace < 9
                    colorSpaceIn = brewermap(nSpace,'Set1');
                elseif nSpace >= 9 && nSpace <=12
                    colorSpaceIn = brewermap(nSpace,'Paired');
                elseif nSpace > 12 && nSpace <=14
                    colorSpaceIn = colorSet(zColSet.rainBow14);
                elseif nSpace == 15 
                    colorSpaceIn = colorSet(zColSet.rainBow15);
                elseif nSpace > 15 && nSpace <=21
                    colorSpaceIn = colorSet(zColSet.rainBow21);
                else
                    colorSpaceIn = jet(nSpace);                    
                end
            end
            
            
            if iscellstr(colorSpaceIn) 
                colorSpace = cell2mat(cellfun(@(x)hexColor(x(2:end)),colorSpaceIn,'uniformoutput',0))';
            elseif ischar(colorSpaceIn)
                cIn = sprintf('colorSpace = %s',colorSpaceIn);
                fprintf('Attempting to use the following color space:\n%s\n',cIn);
                eval(cIn); 
            elseif isnumeric(colorSpaceIn)
                if size(colorSpaceIn,2) ~= 3 
                    error('Colorspace should be array of dx3');
                end
                
                colorSpace = colorSpaceIn(1:min(nSpace,size(colorSpaceIn,1)),:);
            else
                error('Not implemented');                
            end
            
            colorIdx = 1:nSpace;
            if size(colorSpace,1) < nSpace
                fprintf('Insufficient color space: %d < %d\n',size(colorSpace,1),nSpace);
                if opts.wrapColors
                    fprintf('Wrapping space');
                    colorIdx = mod(colorIdx-1,length(colorSpace));
                else
                    error('Add more colors or change to wrapping space');
                end
            end
            
            obj.nameSpace = nameSpace;
            obj.nameSpaceMap = containers.Map(nameSpace,colorIdx);
            obj.colorSpace = colorSpace;
            obj.opts = opts;
        end
        function cm = getColorArray(obj,cStr)
            r = nanvalues(obj.nameSpaceMap,cStr);
            
            er = isnan(r);
            if any(er)    
                if obj.opts.errorOnNew
                    error('Color definition missing');                                 
                end
            
                cm = nan(length(r),3);
            
                cm(~er,:) = obj.colorSpace(r(~er),:);
            
                cm(er,:) = obj.baseColor;
            else
                cm = obj.colorSpace(r,:);
            end
        end
    end
end