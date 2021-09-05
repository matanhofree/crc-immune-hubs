function outarr = concatCell_to_cell(carr)
    isdata = cellfun(@(x)~isempty(x),carr);
    
    carr = carr(isdata);
    sz = sum(cellfun(@(x)numel(x).^iscell(x),carr));
    
    outarr = cell(sz,1);
    cnt = 1;
    for i=1:length(carr)        
%         for j=1:length(carr{i})
%             outarr{cnt} = carr{i}{j};
%             cnt = cnt + 1;
%         end
        zdatCell = carr{i};
        
        if iscell(zdatCell)
            clen = numel(zdatCell);
            outarr(cnt:(cnt+clen-1)) = zdatCell(:);
            cnt = cnt + clen;
        else            
            outarr{cnt} = zdatCell;
            cnt = cnt + 1;
        end
    end
    cnt = cnt - 1;
            
    assert(cnt == sz);
end