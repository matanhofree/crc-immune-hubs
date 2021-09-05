function [elbowVal,elbowIdx] = findElbow(xVal,doPlot,revCdf,verbose)
% Finds elbow in smooth histograms based on data and CDF function (?)
%
    if nargin < 2
        doPlot = 0;
    end
    if nargin < 3
        revCdf = 0;
    end        
    if nargin < 4
        verbose = 1;
    end
    
    [D,N] = size(xVal);    
    if (N>D)
        xVal = xVal';
        [D,N] = size(xVal);    
    end
    
    elbowIdx = nan(1,N); 
    elbowVal = nan(1,N);
    for i = 1:N 
        try
            [xj,Fj] = ecdfQuantile(xVal(:,i));
            if revCdf 
                Fj = 1-Fj;
            end
                
            [elbowVal(i),elbowCdf(i)] = knee_pt(Fj,xj);
            
            [~,elbowIdx(i)] = min(abs(xVal(:,i)-elbowVal(i)));
            
            if verbose
                zrnk = sum(xVal(:,i)>elbowVal(i));
                fprintf('Elbow found at %f (rank: %d)\n',zrnk/size(xVal,1),zrnk);
            end
        catch
            fprintf('Skipped a column in knee point calculation\n');
            continue;
        end

        %% 
        if doPlot == 1
            zf = figure;
            histogram(xVal(:,i));
            addLine(elbowVal(i))
        elseif doPlot == 2 
            zf = figure;
            plot(xj,Fj);
            addLine(elbowVal(i));
        end
    end
end

