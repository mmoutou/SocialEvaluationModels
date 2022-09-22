function xLik = vas2Likert(xVAS, limVAS,limLik)
%vas2Likert discretize (mostly pavlovia) VAS rating xVAS to (default 1-to-6) Likert
%   Default pavlovia goes from 1 to 5 ...

try
    limLik;
catch
    limLik = [1,6];
end
try
    limVAS;
catch
    limVAS = [1,5];  % default Pavlovia VAS as of 2021 outputs from 1 to 5
end
if xVAS < limVAS(1) || xVAS > limVAS(2)
    error('xVAS outside scale');
end

if xVAS == limVAS(1)
    xLik = limLik(1);
elseif xVAS == limVAS(2)
    xLik = limLik(2);
else
    f = (xVAS - limVAS(1))/(limVAS(2)-limVAS(1));
    xL = limLik(1)-1 + f*(limLik(2)-limLik(1)+1);
    xLik = ceil(xL);
end

return;

