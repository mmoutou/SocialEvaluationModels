function natP = tr2nat_SEL_i(P)
% tr2nat_SEL_i - to bring param for ENTIRE SEL data as per selbl03
%                 to native space. P and natP are structures w named fields,
%                 or corresponding vectors, which must be in order
%     1       2         3       4         5       6       7        8     9    10   11      12      
%  posiSelf posiOther dEvSelf dEvOther aEvSelf aEvOther alphaPrec genLR repLR wp0   wAttr  mem    
%__________________________________________________________________________

% Complicated if statement to transform inputted parameters 
if isempty(P)   % if delibarately empty, just provide template w fields in correct order
    natP.posiSelf = nan;   natP.posiOther = nan;  
    
    natP.dEvSelf = nan;    natP.dEvOther = nan;        
    natP.aEvSelf = nan;    natP.aEvOther = nan;   
    
    natP.alphaPrec = nan;  
    
    natP.genLR = nan;     natP.repR = nan;   
    
    natP.wp0 = nan;       natP.wAttr = nan ; 
    natP.mem = nan;   
    
elseif isfield(P,'aEvSelf')    % using this to recognise structure input
    natP =  P;  % allocate space more efficiently and cater for wp0 and wAttr.
    % first, logit-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pTr = P.posiSelf;       natP.posiSelf = 1/(1+exp(-pTr));   
    pTr = P.posiOther;      natP.posiOther = 1/(1+exp(-pTr));       
    pTr = P.genLR;          natP.genLR = 1/(1+exp(-pTr));   
    pTr = P.repLR;          natP.repR = 1/(1+exp(-pTr));   
    pTr = P.mem;            natP.mem = 1/(1+exp(-pTr));   
    % then, log-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    natP.dEvSelf = exp(P.dEvSelf);   
    natP.dEvOther = exp(P.dEvOther);        
    natP.aEvSelf = exp(P.aEvSelf);   
    natP.aEvOther = exp(P.aEvOther);        
    natP.alphaPrec = exp(P.alphaPrec);    
    % scaled-logit transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %     elseif strcmp(field{i},'desBias')
    %         pTr =  P.desBias;  natP.desBias = -1 + 2/(1+exp(-pTr)); 
    % Already copied, already in native space ~~~~~~~~~~~~~~~~~~~~~~~~~~
%     natP.wp0 = P.wp0;  
%     natP.wAttr = P.wAttr;  
else  % 
    if length(P) ~= 12;     error('length(P) ~= 12');    end
    
    natP =  P;  % allocate space more efficiently and cater for wp0 and wAttr.
    %    1       2         3       4         5       6       7        8     9    10   11      12      
    % posiSelf posiOther dEvSelf dEvOther aEvSelf aEvOther alphaPrec genLR repLR wp0 wAttr  mem        
    % first, logit-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    natP(1) = 1/(1+exp(-P(1)));   
    natP(2) = 1/(1+exp(-P(2)));   
    natP(8) = 1/(1+exp(-P(8)));   
    natP(9) = 1/(1+exp(-P(9)));   
    natP(12) = 1/(1+exp(-P(12)));   
    % then, log-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    natP(3:7) = exp(P(3:7));  
    % Already copied, already in native space ~~~~~~~~~~~~~~~~~~~~~~~~~~
%     natP(10:11) = P(10:11);    
end
  
return;
