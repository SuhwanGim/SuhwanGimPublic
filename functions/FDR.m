function [pID,pN] = FDR(p,aFDR)
% FORMAT [pID,pN] = FDR(p,aFDR)
% 
% p     - vector of p-values
% aFDR  - alpha_FDR, aka q, the desired False Discovery Rate level
%
% pID - p-value threshold based on independence or positive dependence
% pN  - Nonparametric p-value threshold
%______________________________________________________________________________
% Author: T. Nichols
% Version: http://github.com/nicholst/matlab/tree/$Format:%h$
%          $Format:%ci$


p = p(isfinite(p));  % Toss NaN's
p = sort(p(:));
V = length(p);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));

pID = p(max(find(p<=I/V*aFDR/cVID)));
if isempty(pID)
  pID=0;
end

pN = p(max(find(p<=I/V*aFDR/cVN)));
if isempty(pN)
  pN=0; 
end
