function [xj,Fj,ecdfScoreOut] = ecdfQuantile(cdfSource)

    cdfSource(isnan(cdfSource)) = [];
    cdfSource(isinf(cdfSource)) = [];
    [Fi,xi] = ecdf(cdfSource);

    if length(xi) > 4

        n = length(xi)-1;
        xj = xi(2:end);
        Fj = (Fi(1:end-1)+Fi(2:end))/2;
        xj = [xj(1)-Fj(1)*(xj(2)-xj(1))/((Fj(2)-Fj(1)));
               xj;
               xj(n)+(1-Fj(n))*((xj(n)-xj(n-1))/(Fj(n)-Fj(n-1)))];
        Fj = [0; Fj; 1];
        
        if nargout > 2
            ecdfScoreOut =  interp1(xj,Fj,cdfSource,'linear','extrap')';
        end
    else
        error('Too few values to derive estimate');
%        ecdfScoreOut = zeros(size(targetScore));
    end


end