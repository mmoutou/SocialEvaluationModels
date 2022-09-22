function AHub = AHubClassLrn(allP, AttrVal)
% AHubClassLrn - A map for probability of positive return.
%
%   allP must include, for use here: wAttr, w0, resNHub, retLevN .
%   AttrVal best reserved for generative process, otherwise omit and use 
%   defaults as per allP.AttrV0

wAttr = allP.wAttr;               w0 = allP.w0; 
resNHub = allP.resNHub;
retLevN = allP.retLevN-1;

Ioff = 0.5;
try 
    if isempty(retLevN); retLevN=2; end
catch
    retLevN = 2;
end
try
    AttrVal;
catch
    AttrVal = allP.AttrV0; % [0.05 0.95]; % for example ...
end
% AHub : Actual positive return depends only on trueAttr
AHub = zeros(retLevN,  resNHub); 

for kA = 1:resNHub
        Attr = AttrVal(kA);
        % Py of LOW return
        AHub(1,kA) = 1- 1/(1+exp(w0+(Attr-Ioff)*wAttr));         
end
AHub(2,:) = 1-(AHub(1,:)); % Py of high (pos) return encoded under high index = 2

return;

