function [midPFair, info] = midpfair(par,res)
%MIDPFAIR provide a grid of 'choice A vs B' probabilities
%   depending on [wH,wS,w0], grid resolution res. 

try wH = par(1);  catch, wH=3; end
try wS = par(2);  catch, wS=5; end
try w0 = par(3);  catch, w0=0; end
try resolN = res; catch, resolN = 4; end
Ioff= 0.5; % just an offset to center the logisitic function below

info = [wH,wS,w0,resolN,Ioff];
midPFair = zeros(resolN,resolN);
midGrid = 1/(2*resolN)+(0:(resolN-1))/resolN;  % grid of midpoints at specified resolution
kH = 1:resolN;
for kS = 1:resolN  
      midPFair(kH,kS) = 1./(1+exp(w0+(midGrid(kH)-Ioff)*wH + (midGrid(kS)-Ioff)*wS)); 
end

return;

