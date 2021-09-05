function outW = weightedNNcontrast(W)

    nD = size(W,2);
    outW = W;
    for i = 1:nD
        zIdx = trueV(i,nD);
        outW(:,i) = W(:,zIdx).*(log1p(W(:,zIdx)) - log1p(nanmax(W(:,~zIdx),[],2)));        
    end
        
end