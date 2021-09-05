function outS = mergeStringPair(strtemp,xv,yv)

    if nargin == 2
        yv = xv;
        xv = strtemp;       
        strtemp = [];
    end
    
    if ischar(xv)
        xv = { xv };
    end
     
    if ischar(yv)
        yv = { yv };
    end
  
    if iscategorical(xv)
        xv = cellstr(xv);
    end
    
    if iscategorical(yv)
        yv = cellstr(yv);
    end
    
    xv = xv(:);
    yv = yv(:);
    
    xN = length(xv);
    yN = length(yv);                
    
    if isempty(strtemp)
        if iscell(xv)
            strtemp = '%s';
        else
            if all(xv == floor(xv))
                strtemp = '%d';
            else
                strtemp = '%f';
            end
            xv = num2cell(xv);           
        end
        
        if iscell(yv)
            strtemp = [ strtemp '_%s'];
        else
            if all(yv == floor(yv))
                strtemp = [ strtemp '_%d' ];
            else
                strtemp = [ strtemp '_%f' ];
            end
            yv = num2cell(yv);
        end
    else
        if ~iscell(xv)
            xv = num2cell(xv);
        end

        if ~iscell(yv)
            yv = num2cell(yv);
        end
    end
    
    if xN ~= yN 
        if xN == 1
            xvE = cell(yN,1);
            xvE(:) = xv;
            
            xv = xvE;
        elseif yN == 1
            yvE = cell(xN,1);
            yvE(:) = yv;        
            
            yv = yvE;
        else             
            error('Dim missmatch');
        end
    end      
    
    outS = cellfun(@(x,y)sprintf(strtemp,x,y),xv,yv,'uniformoutput',0);

end