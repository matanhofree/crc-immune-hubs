function outL = addLine(xPos,yPos,varargin)

    xl = xlim();    
    yl = ylim();
    
    switch nargin
        case 0
            outL = line([ min(xl(1),yl(1)) max(xl(2),yl(2)) ],[ min(xl(1),yl(1))  max(xl(2),yl(2))],'Color','k','lineStyle','--');
        case 1
            if length(xPos) == 1
                xPos = [ xPos xPos ];
            end
            % line([ xPos(1) min(xl(1),yl(1))],[ xPos(2) max(xl(2),yl(2))],'Color','k','lineStyle','--');
            outL = line(xPos,yl,'Color','k','lineStyle','--')
        case 2
            if (isempty(xPos))
                % xPos =  [ min(xl(1),yl(1))  max(xl(2),yl(2)) ];
                xPos = xl;
            end
            
            if (isempty(xPos))
                % yPos =  [ min(xl(1),yl(1))  max(xl(2),yl(2)) ];
                yPos = yl;
            end

            
            if length(xPos) == 1
                xPos = [ xPos xPos ];
            end
            if length(yPos) == 1
                yPos = [ yPos yPos ];
            end
            
            outL =  line(xPos,yPos,'Color','k','lineStyle','--');
            
        otherwise
           error('Too many variables');
    end

end