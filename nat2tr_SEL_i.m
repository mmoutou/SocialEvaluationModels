function transfP8 = nat2tr_SEL_i(par8,vecInOut)
% transfP8 = nat2tr_SEL_i(par8)  : Transform up-to-12 param vector
%            for multiple block Social Evaluation Learning task
%            (Carlisi ... Stothard) 
%            '***8' refers to the (mostly) 8-block experiment
%  Parameters in native space - can be vector or structure:
%         1          2        3      4         5        6       7        8     9    10   11     12    
%      posiSelf  posiOther  dEvSelf dEvOther aEvSelf aEvOther alphaPrec genLR repLR wp0 wAttr   mem
% e.g. [0.75,      0.5,      1,      1,       2,       2,      5,      0.25,  0.8,  0,  6,     0.9998];
try 
    vecInOut; 
catch
    vecInOut = 0;   % if missing, assume 
end
if ~isempty(par8)
    transfP8 = par8;  
    if vecInOut
        transfP8(1)= log(par8(1)/(1-par8(1)));
        transfP8(2)= log(par8(2)/(1-par8(2)));
        transfP8(3:7)= log(par8(3:7));
        transfP8(8)= log(par8(8)/(1-par8(8)));
        transfP8(9)= log(par8(9)/(1-par8(9)));
        %  10 is alread in -INF to INF
        %  11 is alread in -INF to INF
        transfP8(12)= log(par8(12)/(1-par8(12)));
    else
        transfP8.posiSelf  = log(par8.posiSelf/(1-par8.posiSelf));
        transfP8.posiOther = log(par8.posiOther/(1-par8.posiOther));
        transfP8.dEvSelf= log(par8.dEvSelf);
        transfP8.dEvOther= log(par8.dEvOther);
        transfP8.aEvSelf= log(par8.aEvSelf);
        transfP8.aEvOther= log(par8.aEvOther);
        transfP8.alphaPrec= log(par8.alphaPrec);
        transfP8.genLR= log(par8.genLR/(1-par8.genLR));
        transfP8.repLR= log(par8.repLR/(1-par8.repLR));
        % transfP8.wp0 = par8.wp0;                     % done already
        % transfP8.wAttr = par8.wAttr;
        transfP8.mem= log(par8.mem/(1-par8.mem));
    end
else % i.e. if par8 is empty, just return a nan-filled relevant structure.
    if vecInOut
        transfP8 = nan(1,12); 
    else
        transfP8.posiSelf  = nan;
        transfP8.posiOther = nan;
        transfP8.dEvSelf= nan;
        transfP8.dEvOther= nan;
        transfP8.aEvSelf= nan;
        transfP8.aEvOther= nan;
        transfP8.alphaPrec= nan;
        transfP8.genLR= nan;
        transfP8.repLR= nan;
        transfP8.wp0= nan;
        transfP8.wAttr= nan;    
        transfP8.mem= nan;        
    end
end

return;  % -------------------------- eof -------------------------------------------------------------
