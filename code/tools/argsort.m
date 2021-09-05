function [idx,zval] = argsort(varargin)
    [zval,idx] = sort(varargin{:});
    idx(isnan(zval)) = [];
    zval(isnan(zval)) = [];    
end