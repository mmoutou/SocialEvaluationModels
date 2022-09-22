function C = CHubLrn(allP)
% CHubLrn - gives: desirability of positive outcome
% desPos: currency (scaling) of desirability for pos. outcome;
%         allP must include it. 

weights = [-1,1]'  ; % rows unfair, fair (only given at one timestep)

C  =  [ [0 0]' weights]*allP.desPos; 

return;

