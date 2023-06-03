function trP = nat2tr_mdp_L_xii(P)
% nat2tr_mdp_L_xii - to bring param for mdp_L_xii, xiii from native to 
%                  transformed space. natP is a structure w named fields.
%                  If P is deliberately empty, provide a standard structure
%                  with all the fields set to nan.
%__________________________________________________________________________

% Complicated if statement to transform inputted parameters.
if ~isempty(P)
    field = fieldnames(P);  
    for i = 1:length(field)
        % first, log-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if strcmp(field{i},'dInitEv')
            trP.dInitEv = log(P.dInitEv);   
        elseif strcmp(field{i},'aInitEv')
            trP.aInitEv = log(P.aInitEv);   
        elseif strcmp(field{i},'alphaPrec')
            trP.alphaPrec = log(P.alphaPrec);        
    % logit-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        elseif strcmp(field{i},'posi0')
            pNat = P.posi0;       trP.posi0 = log(pNat/(1-pNat)); 
    %     elseif strcmp(field{i},'lrnR')  
    %         pNat = P.lrnR;       trP.lrnR = log(pNat/(1-pNat));  
        elseif strcmp(field{i},'mem')  
            pNat = P.mem;       trP.mem = log(pNat/(1-pNat));  
        % and scaled-logit transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %     elseif strcmp(field{i},'desBias')
    %         pNat =  P.desBias;  trP.desBias = log((1+pNat)/(1-pNat));
        % Already in native space: ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        elseif strcmp(field{i},'wp0') ||  strcmp(field{i},'wAttr') 
           trP.(field{i}) = P.(field{i});  
        else
           % i.e. for w
           % trP.(field{i}) = P.(field{i});  
           error([field{i} ' not catered for here']);
       end
    end
else
    trP.posi0  = nan;
    trP.dInitEv= nan;
    trP.aInitEv= nan;
    trP.alphaPrec= nan;
    trP.wAttr= nan; 
    trP.wp0= nan;   
    trP.mem= nan;     
end
  

return;
