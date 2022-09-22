function dHub = dHubClassLrn(allP)
%DHUBCLASSLRN form Hub d map - priors over Attr, so that the 
%  alpha and beta are each at least 1 and pAttr is the central tendency.
%  allP must include, for use here: pAttr, dInitEv, resNHub, noiseFloor

if length(allP.AttrV0) ~= 2
    error('allP.AttrV0) ~= 2 not catered for yet');
end

attr0 = allP.attr0; 
dInitAttr = allP.dInitEv;
pLeast = allP.noiseFloor;

frAttr = (attr0 - allP.AttrV0(1)) / (allP.AttrV0(2) - allP.AttrV0(1));
if frAttr <= pLeast 
    frAttr = pLeast; 
elseif frAttr >= (1 - pLeast)
    frAttr = 1 - pLeast;
end

dHub{1} = [ (1-frAttr)*(1+dInitAttr); frAttr * (1+dInitAttr)];

return;

