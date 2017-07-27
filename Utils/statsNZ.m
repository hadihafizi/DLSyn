function [ meanNZ, stdNZU, stdNZL, stdNZ ] = statsNZ( wtMat )
%%
%    wtMat - any matrix, really, but a weight matrix is what was in mind 
%
%
% Returns:
%    meanNZ - the mean over each column for all nonzero entries
%    stdNZU/L - the standard deviation (upper and lower) for all nonzeros entries
%
% Description :
%    This finds the first and 2nd moments of the nonzeros values in the columns 
%   of the given matrix
%
%
% Author   : Zach Tosi
%            Indiana University
%
    wtMat(~isfinite(wtMat)) = 0;
    wtMat(isnan(wtMat)) = 0;
    
    adj = wtMat~=0;
    pops = sum(adj);
    meanNZ = sum(wtMat)./pops;
    temp = bsxfun(@minus, wtMat, meanNZ);
    temp = temp .* adj;
    stdNZU = temp .* (temp>0);
    stdNZL = temp .* (temp<0);
    stdNZU = sqrt(sum(stdNZU.^2) ./ (pops-1));
    stdNZL = sqrt(sum(stdNZL.^2) ./ (pops-1));
    stdNZ = sqrt(sum(temp.^2)./( pops-1));
    
end

