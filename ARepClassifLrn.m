function AOne = ARepClassifLrn(allP)
% AClassifRepDict - How well one report (e.g. Attribution) is to be made - written
%                   for separate-attribution-MDP spoke(s).  
% allP must include, for use here: desBias, resNRepo, Ucor, corLevN, resNRepo.

desBias = allP.desBias;    Ucor=allP.Ucor; corLevN=allP.corLevN;
resNRepo = allP.resNRepo;
tiny = 1e-16;    % probabilities of this order are zero for all intends and purposes.
if corLevN == 3
      corGrid = [0.2, 0.5, 0.8];
else
      corGrid = 1/(2*corLevN)+(0:(corLevN-1))/corLevN;
end

%% Calculate AOne : How well Intent is to be reported.
%             (row)          (col)         (page)                  
%         report quality    trueIntent    IntentReport   
AOne =    zeros(corLevN+1,   resNRepo,    resNRepo+1);   % how well HI or SI reported
% At initial report state, all other state components 
% return indifferent/neutral, i.e. probability 1 for corLevN+1 :
AOne(corLevN+1,:,resNRepo+1  ) = 1;    

% Now calc. the resNRepo cols for time 2 on the basis of the trueIntent, 
% e.g. page A3(3, HIlev, SILev, 4, SIrep ) would be high if Intent-lev-to-report=4,
% *shifted by the social desirability bias* ... 
%    To have several levels to account for big or small misattributions. Say we
% call diff of 0 OK, up to 1 lev wrong OKish, the rest just wrong.
for kIntRep = 1:resNRepo         %  Reported Intent level index
     for kIntTru = 1:resNRepo    %  true attribute (as believed!) index
              % first, the unbiased distibution over correctness levels:
              if kIntRep == kIntTru
                  corp = corGrid(corLevN);     % correct case - 'corp' is prob. of 'correct' outcome
              elseif abs(kIntRep - kIntTru) == 1
                  corp = corGrid(corLevN-1);   % approx. OK
              else
                  corp = corGrid(1);           % 'wrong'
              end
              % Now shift it acc. to desirability bias. If desBias is 0, leave corp alone.
              if desBias > tiny         % +ve desir bias is about under-reporting HI, SI
                  if kIntRep >= kIntTru   % Hrep is high enough for pmf to be shifted towards wrong
                     corp = (1-desBias)*corp;
                  else                % Hrep is on the low side, so shift towards correct
                     corp = corp+desBias*(1-corp); 
                  end
              elseif desBias < -tiny  % -ve desir bias is about over-reporting HI, SI.
                  if kIntRep > kIntTru   % Hrep is high enough for pmf to be shifted towards correct
                     corp =  corp-desBias*(1-corp);
                  else                % Hrep is on the low side, so shift towards wrong
                     corp = (1+desBias)*corp; 
                  end
              end  % end if statement that shifts corp according to desirability bias desp.
              
              % Load pmf over correctness levels into AOne 
              %
              AOne(1:corLevN,kIntTru,kIntRep) = noisyBino( corp, Ucor, corLevN ); 
              %
     end  % end loop over all the kIntTrue for this kIntRep
end  % end loop over all kIntRep

return;

