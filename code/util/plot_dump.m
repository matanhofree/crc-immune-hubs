function plot_dump(cfunc,cfig,outPath)

    global printToFigureDump;

    if nargin() < 3
        outPath = 'figures/dump';
    end
    
    if isempty(printToFigureDump) || printToFigureDump == 0
        return;
    end
    
    if ~exist('./figures/dump','dir')
        if printToFigureDump > 1        
            mkdir figures/dump
        else
            message('Figures dump folder not found. Skipping\n');
        end        
    end

    
    if ~exist('cfunc','var') || isempty(cfunc)
        dStack = dbstack;

        if length(dStack) < 1
            cfunc = 'base_plot';
        else
            cfunc = dStack(1:min(3,length(dStack))).name;
            cfunc = strjoin(cfunc);
        end
    end

    try
        if ~exist('cfig','var') || isempty(cfig)
            cfig = gcf();        
        end

        timeStamp = datestr(now,'yyyymmdd_HHMMSS');
        outName = sprintf('%s/%s%s%s',outPath,cfunc,timeStamp) ;   
        print(cfig,outName,'-dpng');        
    catch er
        fprintf('Failed to dump figure\n');
    end
end