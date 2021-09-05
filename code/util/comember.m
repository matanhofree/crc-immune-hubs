function [zia,zib,isCo] = comember(sampleA,sampleB)
    
    uSample = intersect(sampleA,sampleB);
    
    [zia] = ismember(sampleA,uSample);
    [zib] = ismember(sampleB,uSample);
    
    fA = sampleA(zia);
    fB = sampleB(zib);
    
    isCo = isequal(fA(:),fB(:));
%     if ~isCo
%         error('This needs fixing or testing');
%        [~,ziaI,zibI] = intersect(sampleA,sampleB);
%        zia = find(zia);       
%        
%     end
end