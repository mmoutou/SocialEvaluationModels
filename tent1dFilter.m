function [yFilt, wts] = tent1dFilter(inA,xInd,baseL)
%TENT1DFILTER filter a vector inA according to a tent filter
%   length of base in index INTERVALS. So 5 would mean that the very 
%   first tent would go from -2.5 to +2.5 and so would include
%   the first, second and third point of inA:    x-|-x---x---inA(1)---2---3-|-.---...
%                                                  |<--------base 5-------->|
%   The average weight of the filter is 1.
%   return the value at xInd, which should be from 1 to length(inA).

if ~isvector(inA); error('inA should be a 1-D vector/array'); end
sA = size(inA);
Y = inA; if sA(1) ~= 1;  Y = Y'; end  % we'll work with a row vector.
lA = max(sA);
if xInd < 1 || xInd > lA; error('xInd out or range, should be 1...length(inA)'); end

try baseL; catch baseL = ceil(lA/3); end
baseL = baseL + 0.01;   % to avoid such as [1 2 2 1] by ceil etc. below.

% pad out the working vector - easier in case we want to do the whole vector
% etc. in future.
tentMidInd = ceil(baseL/2);
tent = [1:tentMidInd (tentMidInd-1):-1:1];
padN = tentMidInd-1; 
Y    = [zeros(1,padN) Y zeros(1, padN)];
lY   = lA + 2*padN;

% Weight vector initially over the extended working vector:
w = [zeros(1,xInd-1) tent zeros(1,lA-xInd)] ;
% Ensure padding weights are zero:
w(1:padN) = 0;           w((lY-padN+1):end)=0;
% Normalize
wts = w / sum(w); 

yFilt = sum( wts .* Y); 

return;

