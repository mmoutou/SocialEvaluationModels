function AHub = AHubClassLrnPos(allP, AttrVal)
% AHubClassLrnPos - A map for probability of positive return, 
%   much like AHubClassLrn, but based not on w0 but on wp0, and 
%   with the 'high' state resulting in the 'high' outcome, unlike in AHubClassLrn
%   allP must include, for use here: wAttr, w0, resNHub, retLevN,  Ioff .
%   AttrVal best reserved for generative process, otherwise omit and use 
%   defaults as per allP.AttrV0

wAttr = allP.wAttr;               
wp0 = allP.wp0;             
resNHub = allP.resNHub;
retLevN = allP.retLevN-1;
Ioff = allP.Ioff; 

try 
    if isempty(retLevN); retLevN=2; end
catch
    retLevN = 2;
end
try AttrVal; catch, AttrVal = allP.AttrV0; end % [0.05 0.95]; % for example ...

% AHub : Actual positive return depends only on trueAttr
AHub = zeros(retLevN,  resNHub); 

for kA = 1:resNHub
        Attr = AttrVal(kA);
        % Py of LOW return. 
        AHub(1,kA) = 1/(1+exp(wp0+(Attr-Ioff)*wAttr));         
end
AHub(2,:) = 1-(AHub(1,:)); % Py of high (pos) return encoded under high index = 2

return;

