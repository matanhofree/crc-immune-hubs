function outarr = concatCell(carr)
    isdata = cellfun(@(x)~isempty(x),carr);
    
    carr = carr(isdata);
    sz = sum(cellfun(@(x)length(x),carr));
    
    outarr = zeros(sz,1);
    cnt = 1;
    for i=1:length(carr)
        for j=1:length(carr{i})
            try 
                outarr(cnt) = carr{i}(j);
                cnt = cnt + 1;
                continue;
            catch
            end

            try 
                outarr(cnt) = carr{i}{j}(:);
                cnt = cnt + 1;
                continue;
            catch
            end
        end
    end
            
        
end