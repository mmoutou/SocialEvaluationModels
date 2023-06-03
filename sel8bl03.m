function [fitM, Par8, MDP8] = sel8bl03( pTr8, inp8, resp8, modStruc8, repoFn, details)
% sel8bl puts together all 8 blocks of an experiment. pTr8 are the participant param (transf. -> -INF to INF),
%        inp8 the inputs (stimuli) that the participant saw over all 8 blocks, resp8 the responses 
%        they gave. If resp8 (or both resp8 and inp8) are not given, then the function runs in
%        generative mode and simulates an experiment. parPri are the priors on the parameters,
%        which are optional.
%            The function outputs are as per learnReportXX, but now for 8 blocks, plus fitM
%        are the fit measures sum-log-likelihood and a posteriori density
testCode = 0;
try repoFn;  catch, repoFn = 1; end    % repoFn 0 is the AcIn report MDPs, repoFn 1 based on
                                       % full pmf sharpening (softmax-like) ... 
try details; catch, details=1;  end

%% About parameters ------------------------------------------------------------------------------- 
%            **** First, parameters which we need to prioritise  *******
% attrSelf, attrOther  : Initial belief of how where between 'negative' and 'positive' the rater 
%            is likely to be. So 1/3 means 2-to-1 'negative', 2/3 means 2-to-1 'positive'
%            It is necessary to explore / fit this!
% dEvSelf/Other (self / other, possibly) : measure of confidence in attrSelf/Other. 
% aEvSelf/Other, possibly) : strength of belief re. how positive 'positive raters' and 
%             'negative raters' are. A lesser priority to explore, first set aEvSelf=aEvOther.
% alphaPrec    : Decision noise parameter
% genLR    : General between-blocks learning parameter
% repLR    : between-blocks learning param for 'repeat' blocks, eg. between 1st and 2nd 'Mary rates Me'.
% w0       : Overall positivity bias 
%      Next,           **** parameters that can be left alone in the first instance  ****
% wAttr    : Slope between attributes and behaviours. Initially fixed at 6, so that 
%            if the rater is '95% bad', and w0=0, then they will give a good rating only 
%            1/(1+exp(w0+(0.95-0.5)*6 )) = 0.063 of the time.
% lrnR     : Active-inference learning rate within blocks. Leave it alone at 1.
% [desBias  : social desirability parameter - dropped, as can be absorbed in w0 in this experiment.]

% Transform parameters back to native space, in order to work with them for learning etc.
% REM params :      1       2        3      4         5        6       7        8     9    10   11      12
% at most are: posiSelf posiOther dEvSelf dEvOther aEvSelf aEvOther alphaPrec genLR repLR wp0   wAttr   mem  
par8 = pTr8;  
par8(1)=  1/(1+exp(-pTr8(1)));  
par8(2)=  1/(1+exp(-pTr8(2)));  
par8(3:7)= exp(par8(3:7));
genLR = 1/(1+exp(-pTr8(8)));            par8(8) =  genLR;
repLR = 1/(1+exp(-pTr8(9)));            par8(9) =  repLR; 
par8(12)= 1/(1+exp(-pTr8(12))); 

Par8.forFit = pTr8;        

%% About stimuli and conditions --------------------------------------------------------------------------- 
try inp8;  catch, inp8 = [];        end
try resp8; catch, resp8 = [];       end

totBlocks = length(inp8);        
pPosGen=[0.2,0.5,0.8]; try pPosGen=inp8{1}.allCodes{5}; catch ; end   % for neg, neutr, pos rater blocks.
resNRepo =  modStruc8.allP.resNRepo; 
midGrid = 1/(2*resNRepo)+(0:(resNRepo-1))/resNRepo;

MDP8 = {};  Resp8={}; Inp8={};     % These will store all blocks for the new estimation ...
% Cater for generative mode:
if isempty(resp8); resp8={};  for blN=1:8; resp8{blN} = []; end; end

% Numerical condition codes. These can code for any ratee-valence-repetition of a block:
% 110:  selfNeg    120:selfNeu     130: selfPos [rater NOT repeated]
% 111:  selfNeg1   121:selfNeu1    131: selfPos1 [rater WILL be repeated]
% 120:  selfNeg2   122:selfNeu2    132: selfPos2 [rater REPETITION ]
% 210:  otherNeg   220:otherNeu    230: otherPos [rater repeated]
% 211:  otherNeg1  221:otherNeu1   231: otherPos1 [rater be repeated]
% 220:  otherNeg2  222:otherNeu2   232: otherPos2 [rater REPETITION ]

% In Isabel's scheme, there are two 'same reater but Pos' repeats, one of 
% which has to follow the corresp. Neg rater, according to Isabel's sequencing, 
% So  a '4', let's say 'Tom -ve about Mo' comes shortly before a '6', 'Tom  +ve about Mo' 
% and a '1'            'Zara -ve about you' shotly before  a     '3', 'Zara +ve about you',  E.g.: 
%                       Zara-you  Tom-Mo Zara+you Tom+Mo
% These are entries in the form of Isabel's data files:
% inp8{1}.allCodes = {condCode,rater_name,ratee_name,ratee_type,pPosGen};

        
%% Core block of code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
selfStarted = 0;           otherStarted = 0;     % to flag first trial where Self or Other are rated.

% REM }     1          2        3      4         5        6       7         8     9    10   11   12  
% pTr8}  posiSelf  posiOther  dEvSelf dEvOther aEvSelf aEvOther alphaPrec genLR repLR  wp0 wAttr mem 
% whereas selp must be;
%        [posi0,              dInitEv,          aInitEv,        alphaPrec,      ** wAttr,w0,**  mem]
pSelfIndices = [1,3,5,  7,11,10,12];           % indices of pTr8 that correspond to entries of rateeP
                                               % for 'Self rated' blocks

DEnd = nan(2,totBlocks);   % posterior beliefs FOR EACH BLOCK
totRatN = max(inp8{1}.allCodes{6}(1,:));   % allCodes{6} should be raterIDs
RaterP   = zeros(1,totRatN);      % coarse measure of overall posterior FOR EACH RATER
blRaterID = inp8{1}.allCodes{6};  % The sequence of raters (with reference to this experiment) for each block.
noiseFl = modStruc8.allP.noiseFloor; 
pHi = 1-noiseFl;
pLo = noiseFl;

for blN = 1:totBlocks
    cond = inp8{blN}.condCode;
    ratee = blRaterID(2, blN);
    if ratee == 1
       rateeP = par8(pSelfIndices); 
    else
       rateeP = par8([2,4,6, pSelfIndices(4:end)]);
    end

    if blN == 1  % i.e. use the input parameters straight - no learning to do! ====================
        % adjust parameters depending if the Self or Other is rated:
        Par8.trialP = nan(totBlocks,length(rateeP));
        
        Par8.trialP(blN,:) = rateeP;
        selfP = spm_unvec(rateeP,nat2tr_mdp_L_xii([])); % was modStruc8.indexP);
        transfP = nat2tr_mdp_L_xii(selfP);

        MDP8{blN} = spm_mdp_L_xiii(transfP,modStruc8,inp8{blN},resp8{blN},repoFn,details); %main ...
                    % structure of MDPs likelihood values, POMDPs ...  
        if testCode; disp(inp8{1}.allCodes{2});  disp(inp8{1}.allCodes{1});  end             
        
    else  %% ======================= if we have learning to do! ===================================
        % Fill in components from end of previous block iteration.
        dend = MDP8{blN-1}.hubHMM(end).d{1};  % fill in posteriors array.
        DEnd(:,blN-1) = dend / sum(dend);
        prevID = blRaterID(1,blN-1);  % this is pt specific - see toySEL8bl01.m for how it's done...
        thisID = blRaterID(1,blN);    % ... above for prev. block, now for this block.
        prevCond = inp8{blN-1}.condCode;
        % How many repetitions of that rater have we seen:
        raterReps = max(1, rem(prevCond,10));
        % Accumulate the average, approximate beliefs re. this rater:
        RaterP(prevID) = (RaterP(prevID)*(raterReps-1) + DEnd(1,blN-1))/raterReps;
  
        % Raters that are not about to be repeated, seen so far.
        % First, all up to the max so far have been seen, and in seq:
        ratSoFar = 1:max(blRaterID(1,1:blN));
        % Now separate the one about to be encountered, if present:
        restRat = ratSoFar(ratSoFar(1,:) ~= thisID);  
        
        % Calculate updates for posiSelf and posiOther based on the unrepeated raters --------------
        % First, make space & initialise w posiSelf and posiOther. The first col. is for
        % the initial conditions, the last is used if the rater thisID has been repeated.
        cnt = 1; 
        posiUpd = nan(2,length(restRat)+2); 
        posiUpd(:,cnt) = par8(1:2); % rows posiSelf, posiOther.
        dEviUpd =  posiUpd;  % rows dEvSelf, dEvOther.  LEAVE
                             % aEvSelf, aEvOther ALONE FOR NOW.
        dEviUpd(:,cnt) = par8(3:4);                  
      
        for rater = restRat
            
            ratee = blRaterID(2, find(blRaterID(1,:)==rater,1));
            
            % prediction error will be used to correct both posiSelf and posiOther.
            % Note that each rater should have only one category of ratee, Self or Other, so
            % PE = Posterior-for-rater - prior-for-corresp.-Self/Oth
            upd = genLR * (RaterP(rater) - posiUpd( ratee, cnt)); 
            
            % Now update!
            cnt = cnt + 1;
            posiUpd(:,cnt) =  posiUpd(:, cnt-1) + upd ;
            tooHi = posiUpd(:,cnt) > pHi ;                    posiUpd(tooHi,cnt) = pHi; 
            tooLo = posiUpd(:,cnt) < pLo ;                    posiUpd(tooLo,cnt) = pLo; 
            
            % The following is embarassingly rough ...
            dEviUpd(:,cnt) = dEviUpd(:,cnt-1);
            dEviUpd(ratee,cnt) = dEviUpd(ratee,cnt-1) + genLR ; 
                       
        end
        
       % Updates for posiSelf and posiOther if rater has been repeated -----------------------
       % Check if rater has been repeated, and if so when last:
        raterLastSeen = find(blRaterID(1,1:(blN-1))==thisID, 1, 'last' );
        if ~isempty(raterLastSeen)
            ratee = blRaterID(2, raterLastSeen); 
            posiUpd(:,end) = posiUpd(:,end-1);
            posiUpd(ratee,end) = repLR * RaterP(thisID) + (1-repLR)* posiUpd(ratee,end);
            
            dEviUpd(:,end) = dEviUpd(:,end-1);
            dEviUpd(ratee,end) = dEviUpd(ratee,end) + repLR;
            
        % =================== Now to update parameter vector! ==================================    
        % REM }     1          2        3      4         5        6       7     8     9    10   11   12     
        % par8}  posiSelf  posiOther  dEvSelf dEvOther aEvSelf aEvOther alphaPrec genLR repLR  wp0 wAttr lrnR 
        % whereas selp must be;
        %        [attr0,              dInitEv,         aInitEv,         alphaPrec,      ** wAttr,w0,**  lrnR ]
         
            rateeP(1) = posiUpd(ratee,end);
            rateeP(2) = dEviUpd(ratee,end);
        else  % if the rater-ratee have not been repeated, then the end col of posiUpd and dEviUpd
              % should be nan and we should use entries from the previous one:
            ratee = blRaterID(2, blN);   
            rateeP(1) = posiUpd(ratee,end-1);
            rateeP(2) = dEviUpd(ratee,end-1);
        end
        if testCode
           disp('posiSelf/Oth updates:'); disp(posiUpd); 
           disp('d Evidence Self/Oth updates:'); disp(dEviUpd); 
        end
        
        Par8.trialP(blN,:) = rateeP;
        selfP = spm_unvec(rateeP,nat2tr_mdp_L_xii([])); % was modStruc8.indexP);
        transfP = nat2tr_mdp_L_xii(selfP);
        
        MDP8{blN} = spm_mdp_L_xiii(transfP,modStruc8,inp8{blN},resp8{blN},repoFn,details); %main structure 
             % for MDPs, including likelihood values, POMDPs ...

    % Rationale for learning attrSelf, attrOther between blocks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % 1. Before calculating likelihood, consider all the raters that have been 
    %    seen so far, and apply learning anew, based on the initial attrSelf and attrOther,
    %    f0s and f0o for short. The learning rate genLR is appliced per individual. 
    % 2. We assume that the posteriors about each rater fully capture the attitude of that
    %    rater for that block; that the respective priors have been washed out by the 20 trials.
    % 3. So, if we have raters u z1 y1 x y2 z2 v w, then just before y2, say self is being rated, 
    %    apply all the prediction errors that would have been encountered as follows:
    %    f0s' <- f0s'+genLR*D(Pk,f0*);   for k in u, z, x,  i.e., NOT y1   and where f0* 
    %    is the f0 appropritate to that D. So we need to update f0s AND f0o for u, z, x
    % 4. when we are about to calculate y2, we do
    %    f0sTemp <- (1-repLR)*f0s + repLR*Py1 
    % 5. When about to do z2, 3. is applied to u,y (now (Py1+Py2)/2 ), x, then 4. from z1.
    %    (z1 is by now some way in the past but it doesn't matter).
    % 6. For v, 3. is applied to u, z, mean y, mean x, in that order.
    % Note that the prediction errors are calculated ANEW with respect to f0s, e.g.
    % f0s has been updated when D((Py1+Py2)/2,f0s) is calculated, but the posteriors are not,
    % following the approximation they are independent of the priors that were used to derive them,
    % ie. that the 20 trials have made minor changes in the estimate of the priors used to derive Py2 etc.
    % immaterial.      
    
        
     % just for completeness and checking - not needed for learning from block to block
       if testCode
        disp(['rateeP updated by between-block learning after block ' num2str(blN)]); 
        disp(rateeP); 
 
        if blN == totBlocks   
         dend = MDP8{blN}.hubHMM(end).d{1};  % fill in posteriors array.
         DEnd(:,blN) = dend / sum(dend);
         ID = inp8{blN}.raterID;    % this is pt specific - see toySEL8bl01.m for how it's done.
         % How many repetitions of that rater have we seen:
         raterReps = max(1, rem(cond,10));
         % Accumulate the average, approximate beliefs re. this rater:
         RaterP(ID) = (RaterP(ID)*(raterReps-1) + DEnd(1,blN))/raterReps;    
        end  % blN == totBlocks  % just for completeness and checking
       end 
    end
    
end

if ~isempty(modStruc8.priPar) 
% REM e.g.:           1       2        3      4         5        6       7        8     9    10   11      12      
% fields:        posiSelf posiOther dEvSelf dEvOther aEvSelf aEvOther alphaPrec genLR repLR wp0   wAttr  lrnR    
% modStruc8.priPar=[[1.01,  1.01,     1.01,    1.01,  1.01,    1.01,     1.01,   1.01, 1.01, 10,    1.2,  1.01]; ...  % A
%                     [1.01,  1.01,      2,       2,      2,       2,        2,    1.01, 1.01, 10,    5.8,  1.01]; ...  % B
%                     [ 0,     0,        0,       0,      0,       0,        0,     0,   0,   -46,     0,    0   ]; ...  % lo
%                     [ 1,     1,        100,    100,    100,     100,     100,     1,   1,    46     100,   1   ]];     % hi
%                                    <-  max             at 1                >            <SD ~10> <max at 4> 
   priPar = modStruc8.priPar;
   logPri = log( dbetasc(par8,priPar(1,:),priPar(2,:),priPar(3,:), priPar(4,:)));     
else
   logPri = 0;
end

if details == 0   % then just calculate the posterior density, or sum LL if 
                  % no prior over parameters has been specified:
  fitM = logPri; 
  for blN=1:totBlocks; fitM = fitM + MDP8{blN}.sum; end
else  % record LL for each block, sumLL, and sum ln posterior density:
  fitM.LL = zeros(1,totBlocks);
  for blN=1:totBlocks; fitM.LL(blN) = MDP8{blN}.sum; end
  fitM.sLL = sum(fitM.LL);
  fitM.logPri = logPri;
  fitM.pDens = sum(logPri) + fitM.sLL; 
end


return; 

%% =================================  eof =====================================
