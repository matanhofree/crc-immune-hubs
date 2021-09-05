function structOut = structpack(structInput,varargin)
% Pack a structure into component variables
    if (~isempty(structInput))
        if (isstruct(structInput))
            structOut = structInput;
        else
            zname = inputname(1);
            structOut = struct();
            structOut.(zname) = structInput;
        end
    else
        structOut = struct();
    end
    
    for i = 1:length(varargin)
        zname = inputname(i+1);
        if (~isempty(zname))
            structOut.(zname) = varargin{i};
        end
    end
        
end