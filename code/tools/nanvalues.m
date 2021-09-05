function [res,isKeyMissing,missing] = nanvalues(map,array)
    
    isKeyMissing = false(size(array));
    missing = {};

    cnt_miss = 1;
    
    if (ischar(array)) % Treat single chars linge cell array of chars
        array = { array };
    end
    
    if (iscell(array))                        
        if (any(strcmp(map.ValueType,{'char' 'any'})))       
            res = cell(size(array));
            for i = 1:numel(array)
                if (map.isKey(array{i}))
                    res{i} = map(array{i});
                else
                    isKeyMissing(i) = 1;
                    missing{cnt_miss} = array{i};
                    cnt_miss = cnt_miss + 1;
                end
            end
        else
            res = nan(size(array));
            for i = 1:numel(array)
                if (map.isKey(array{i}))
                    res(i) = map(array{i});
                else
                    isKeyMissing(i) = 1;
                    missing{cnt_miss} = array{i};
                    cnt_miss = cnt_miss + 1;
                end
            end
        end
    else
        
        cIdx = 1;
        try      
            if ~(strcmp(map.ValueType,'char') == 1 || strcmp(map.ValueType,'any') == 1)
                res = nan(size(array));
                elemCnt = numel(array);
                
                for i = 1:elemCnt
                    if (map.isKey(array(i)))
                        cRes = map(array(i));
                        res(i) = cRes;
                        % res(cIdx:(cIdx+length(cRes)-1),:) = cRes';
                        % res(cIdx:(cIdx+length(cRes)-1)) = cRes';
                        % cIdx = cIdx + length(cRes);
                    else
                        isKeyMissing(i) = 1;
                        missing{cnt_miss} = array(i);
                        cnt_miss = cnt_miss + 1;
                    end
                end
            else               
                res = cell(size(array));
                               
                for i = 1:numel(array)
                    if (map.isKey(array(i)))
                        res{i} = map(array(i));
                    else
                        isKeyMissing(i) = 1;
                        missing{cnt_miss} = array(i);
                        cnt_miss = cnt_miss + 1;
                    end
                end
                
            end
        catch e
            error(1,'Error: failed attempting to use cell output format\n');
   
        end         
    end
    
end