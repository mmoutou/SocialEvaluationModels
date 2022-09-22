function natP = tr2nat_mdp_L_xi(P)
% tr2nat_mdp_L_xi - to bring param for mdp_L_xi to iv from transformed to 
%                  native space. P and natP are structures w named fields.
% test/demo with tr2nat_mdp_L_xi(nat2tr_mdp_L_xi(nativePar))
%__________________________________________________________________________

% Complicated if statement to transform inputted parameters 
field = fieldnames(P);  
for i = 1:length(field)
    % first, log-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if strcmp(field{i},'dInitEv')
        natP.dInitEv = exp(P.dInitEv);   
    elseif strcmp(field{i},'aInitEv')
        natP.aInitEv = exp(P.aInitEv);        
    elseif strcmp(field{i},'alphaPrec')
        natP.alphaPrec = exp(P.alphaPrec);    
    % probably won't be used, but harmless:
    elseif strcmp(field{i},'Ucor')
        natP.Ucor = exp(P.Ucor);  
    % logit-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif strcmp(field{i},'attr0')
        pNat = P.attr0;       natP.attr0 = 1/(1+exp(-pNat));  
    elseif strcmp(field{i},'lrnR')  
        pNat = P.lrnR;      natP.lrnR = 1/(1+exp(-pNat));   
    elseif strcmp(field{i},'mem')  
        pNat = P.mem;       natP.mem = 1/(1+exp(-pNat));   
    % and scaled-logit transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif strcmp(field{i},'desBias')
        pNat =  P.desBias;  natP.desBias = -1 + 2/(1+exp(-pNat)); 
    % Assume everything else Untransformed ~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
       natP.(field{i}) = P.(field{i});  
   end
end
  

return;
