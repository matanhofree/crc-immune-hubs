function varargout = minmax(a, type, whis)
%MINMAX Returns minimum and maximum value in the given array
%
% [minval maxval] = minmax(a)
% lims = minmax(a);
% lims = minmax(a, type);
% lims = minmax(a, type, w);
%
% Computes the minimum and maximum value in entire array (all dimensions).
%
% Input variables:
%
%   a:      numeric array
%
%   type:   'all':          absolute minimum and maximum (default)   
%           'noout':        discards outliers
%           'center':       centers on zero
%           'centernoout':  centers on 0 and eliminates outliers
%           'expand':       wider that 'all' by a specifed fraction
%
%   w:      for no-outlier version, factor defining an outlier.  Point is
%           considered an outlier if larger than q3+w*(q3-q1) or smaller
%           than q1-w*(q3-q1) [1.5]  
%           
%           for expand version, fraction of actual range to add onto each
%           end [0.1]
%
% Output variables:
%
%   minval: minimum value in a
%
%   maxval: maximum value in a

% Copywrite 2005 Kelly Kearney

if nargin == 1
    type = 'all';
end

switch type
    case 'all'
        minval = min(a(:));
        maxval = max(a(:));
    case 'noout'
        if nargin < 3
            whis = 1.5;
        end
        pctiles = prctile(a(:),[25 75]);
        q1 = pctiles(1);
        q3 = pctiles(2);
        vhi = q3+whis*(q3-q1);
        vlo = q1-whis*(q3-q1);
        minval = min(a(a > vlo));
        maxval = max(a(a < vhi));
        if isempty(minval)
            minval = vlo;
        end
        if isempty(maxval)
            maxval = vhi;
        end
    case 'center'
        temp = max(abs([min(a(:)) max(a(:))]));
        minval = -temp;
        maxval = temp;
    case 'centernoout'
        if nargin < 3
            whis = 1.5;
        end
        pctiles = prctile(a(:),[25 75]);
        q1 = pctiles(1);
        q3 = pctiles(2);
        vhi = q3+whis*(q3-q1);
        vlo = q1-whis*(q3-q1);
        temp = max(abs([min(a(a > vlo)) max(a(a < vhi))]));
        minval = -temp;
        maxval = temp;
    case 'expand'
        if nargin < 3
            whis = 0.1;
        end
        minval = min(a(:));
        maxval = max(a(:));
        da = maxval - minval;
        minval = minval - whis*da;
        maxval = maxval + whis*da;
    otherwise
        error('Unrecognized option: %s', type);
end
        
       
if nargout == 2
    varargout{1} = minval;
    varargout{2} = maxval;
elseif nargout == 1
    varargout{1} = [minval maxval];
elseif nargout == 0
    [minval maxval]
else
    error('Wrong number of output arguments');
end