function trP = nat2tr_mdp_L_xi(P)
% nat2tr_mdp_L_xi - to bring param for mdp_L_xi from native to 
%                  transformed space. natP is a structure w named fields.
%__________________________________________________________________________

% Complicated if statement to transform inputted parameters.
field = fieldnames(P);  
for i = 1:length(field)
    % first, log-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if strcmp(field{i},'dInitEv')
        trP.dInitEv = log(P.dInitEv);   
    elseif strcmp(field{i},'aInitEv')
        trP.aInitEv = log(P.aInitEv);   
    elseif strcmp(field{i},'alphaPrec')
        trP.alphaPrec = log(P.alphaPrec);     
    % probably won't be used, but harmless:
    elseif strcmp(field{i},'Ucor')
        trP.Ucor = log(P.Ucor);      
% logit-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif strcmp(field{i},'attr0')
        pNat = P.attr0;       trP.attr0 = log(pNat/(1-pNat)); 
    elseif strcmp(field{i},'lrnR')  
        pNat = P.lrnR;       trP.lrnR = log(pNat/(1-pNat));  
    elseif strcmp(field{i},'mem')  
        pNat = P.mem;       trP.mem = log(pNat/(1-pNat));  
    % and scaled-logit transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif strcmp(field{i},'desBias')
        pNat =  P.desBias;  trP.desBias = log((1+pNat)/(1-pNat));
    % Assume every thing else Untransformed ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
       trP.(field{i}) = P.(field{i});  
   end
end
  

return;
