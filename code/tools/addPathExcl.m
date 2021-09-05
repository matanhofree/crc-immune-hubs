function oldpath = addPathExcl(inPath,verbose)
    
    if ispc()
        inPath = addpathPlus(inPath);
        inPath = inPath{:};
    end

    if isunix() & ~ismac()
         %%
        [ ~, prepPath] = system(sprintf('find %s -type d',inPath));
        prepPath = textscan(prepPath,'%s','delimiter','\n');
        prepPath = prepPath{1};
    elseif ispc()
        prepPath = genpath(inPath);
        prepPath = textscan(prepPath,'%s','delimiter',';');
        prepPath = prepPath{1}; 
    elseif ismac()
        
        [ cout, prepPath] = system(sprintf('find %s -type d',inPath));
        
        if cout ~= 0 || isempty(prepPath)   
            [ cout, prepPath] = system(sprintf('/usr/local/bin/gfind %s/ -type d',inPath));
        end
        
        if cout ~= 0 || isempty(prepPath)            
            prepPath = genpath(inPath);
            prepPath = textscan(prepPath,'%s','delimiter',';');
            prepPath = prepPath{1}; 
        else
            prepPath = textscan(prepPath,'%s','delimiter','\n');
            prepPath = prepPath{1};

        end        
    else
        prepPath = genpath(inPath);
        prepPath = textscan(prepPath,'%s','delimiter',':');
        prepPath = prepPath{1};
    end
    
    %% Drop hidden dirs (e.g. git repositories)    
    if ispc()
        dropRepoDirs = strgrep(prepPath,'.*\\\..+');
        dropRepoDirs = dropRepoDirs | strgrep(prepPath,'.*\/\..+');
    else
        dropRepoDirs = strgrep(prepPath,'.*\/\..+');    
    end
    prepPath(dropRepoDirs) = [];
    
    %% Drop dir without m files    
    dropNonM = cellfun(@(xdir)isempty(dir([ xdir filesep '*.m'])),prepPath);
        %& ...
        % cellfun(@(xdir)isempty(dir([ xdir filesep '@*'])),prepPath);
    prepPath(dropNonM) = [];
    
    %% Drop class dirs    
    if ispc()
        dropRepoDirs = strgrep(prepPath,'.*\\@.+');        
    else
        dropRepoDirs = strgrep(prepPath,'.*\/@.+');    
    end
    prepPath(dropRepoDirs) = [];
    
%%        
    if exist('verbose','var') && verbose
        cellfun(@(x)disp(x),prepPath);
    end
%%
    if ispc()
        prepPath(:,2) = {';'};
    else
        prepPath(:,2) = {':'};    
    end
    
    prepPath = prepPath';
    prepPath = [ prepPath{:} ];
    oldpath = addpath(prepPath);
    
end