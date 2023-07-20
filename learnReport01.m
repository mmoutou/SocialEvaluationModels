function [MDPs, modelStructure, Inp, Resp] = learnReport01(selfp, pRet, toPlot)
% learnReport01 - reporting chance of positive outcome, updating
%        beliefs based on actual return, AND attribution reporting.
%        Version to use good-bad classification plus A learning, after
%        Pike 'exponentially reducing learning rate'.
%        01 has 1 independent reporting MDPs, pRew reporting
%        MDPs draw inputs from the learning-MDP, (up to 06 had 1 attribution MDP).
%        This factorises the solution and means that there is no feedback from
%        'Attr  (attributions) I will/have declared' to 'underlying beliefs about Attr'.
%        toPlot =0: no plots; 1: behaviour & key beliefs; 2: more beh. detail. 3: key neural simulations.
%        6 levels of reporting, 'boring' desBias (desirability bias) could be used 
%        in AClassifAttrLrn), but added intEvRat
% Demo: 
% close all; pRet=0.6; selfp= [0.5,   1,      2,     1/5,   6,   0,   1,      0]; [MDPs,modStruc,Inp,Resp] = learnReport01(selfp,pRet,1);
%                             [attr0,dInitEv,aInitEv,uPrec,wAttr,w0,lrnR, desBias]
% Uses hidden states causing (factorised) outcomes. Here the factorisation is
% explicit, enabling us to model multiple modalities (outcome factors) and
% distinct hidden causes of observation (hidden state factors like what and
% where).
%______________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging, Karl Friston et al

testing = 0; % 0 -> use time-based random gen ...
seed = 661;    % Only used if testing~=0

% set up and preliminaries
%==========================================================================
if ~testing
    rng('shuffle');
else
    rng(seed, 'twister');
    warning('learnReport01 running with testing/debug settings incl rng(seed,''twister''');
end % for reproducible sequences of random numbers; 'shuffle' bases on time.

%% Menu - like items -----------------------------------------------------
% Number of trials (maybe in future per each Other encountered)
trialN = 30; % I think Kate's task has 30 trials per rater, x 6 rater-rated combos?
try selfp;  catch, selfp  = []; end
try pRet; catch, pRet = []; end
if isempty(pRet);  pRet = 0.6;   end  % True positive-return probability

% Use the pRet above to construct 'true' binary class:
trueAttr = 1+1*(pRet(1) > 0.5);  % call these, a little arbitrarily, 'high' and 'low'
% True 'alternatives'; or pseudo-alternatives, as true state is fixed to one of them:
if trueAttr == 1; AttrVal = [pRet 0.75];  else AttrVal = [0.25 pRet];  end

%                               1     2         3     4      5    6      7    8           
%                             [attr0,dInitEv,aInitEv,uPrec,wAttr, w0,  lrnR, desBias]
if isempty(selfp);   selfp  = [0.1,   2,       2,    0.1,    6,   0,   1,    0.0];   end
                                 % REM Joe's reversal paper had attr0,pS0,uH,uS, wAttr,wS,w0,uPi,eta_dg
                                 % uPrec (cf. uPol) is 1/alphaPrec, the expected prior prec.
attr0 = selfp(1);      
dInitEv = selfp(2);   aInitEv = selfp(3);  % give error if these fail :)
try
    alphaPrec = 1/selfp(4); % uPrec (cf. uPol) is 1/alphaPrec, the expected prior prec.
    wAttr=selfp(5); w0 = selfp(6); lrnR= selfp(7);
    %    bias from -1 to +1, to shift pmf; ...]
    desBias = selfp(8);
catch
    alphaPrec = 10; wAttr=5; wS=5; w0=0.5; lrnR=1; desCorr = 1;
end

Ucor = 0.2; % Noise param for noisyBino over wrong ... correct. 0.2 -> care, 2 -> don't care
            % See A3 and A4AttrLrn - How well Attr, SI is to be reported.
Upop = 1.0; % see dClassifLrn -- param for belief of SD of population, reflected
            %     in bluntness of d. So set to 1 for representative mean, 0.2 for representative mode (?).
desPos=1;  % scaling of desirability to get pos return.
desCorr=1;  % scaling of desirability to report correctly, given one's beliefs.

% toPlot = 4 to display neuronal phase precesion / theta-gamma etc.
%       or 1-3 for more basic plots.
try toPlot; catch, toPlot =0; end
% neuroFami = 1 to display familiarity / MMN / other_intent learning - NOT READY.
neuroFamil =0; %try, neuroFamil; catch, neuroFamil =0; end

%% Model specifiation: dimensions of matrices and descriptions of said dimensions/factors
Tsteps1 = 2;        % Each stage of learning-MDP, hmmHub will have 2 Tsteps.
                    %  initial state -> receive-return state.
Tsteps2 = 2;        % Each stage of reporting MDP, will only have just 2 Tsteps.
                    %  initial state -> choose level-to-report
stateFNHub = 1;     % Number of state factors for hub HMM == learning-MDP: trueAttr
stateFNRep = 2;    % N. of state factors for reporting-MDP: beliefState, reportState
% N. of (pseudo) observations of whether we reported well or not. Here using 3 levels
% of correctness of reporting, e.g. quite correct=3, not too bad, ...
% ...at most 1 off in each dimension =2, wrong ie the rest =1
corLevN = 3;   %  corLevN+1 is the 'indifferent/neutral'

resNHub  = 2;   % resoln. level for classifying partner, starting w. e.g.
                % neg (0.1 return) mid (0.45) pos (0.8)
resNRepo = 6;   % Likert-like resolution for reporting positivity or attributions - resolN intervals, w.
                % boundaries (1:(resolN-1))/resolN and midpoints 1/(2*resolN)+(0:(resolN-1))/resolN


midGridRepo = 1/(2*resNRepo)+(0:(resNRepo-1))/resNRepo;  % grid of midpoints at specified resolution
ubGridRepo  = (1:resNRepo)/resNRepo;                   % grid of upper bounds of same.
resolRepo   = 1/resNRepo;

% Generative model parameters to implement Other's policy:
desBias = 0;  % Usually don't use this but w0 to account for overall biases if poss.
              % Can be used (is coded into) in AClassifAttrLrn
% midPFair is a resoln x resoln grid of 'choice A vs B' probabilities depending
% on [wAttr,0,w0]
midPPos = midpfair([wAttr,0,w0],resNHub);

Ioff = 0.5;  % Always use the same offset for the logistic below, for centering
corePFair = 1/(1+exp(w0+(AttrVal(trueAttr)-Ioff)*wAttr ));
[~, corePBin] = min(abs(midGridRepo - corePFair));  % Bin where the intended generative
                                                    % pos-share proportion resides.

% Make a grid of index probabilities to aid noisify A2 (attributing) map. Each
% col of pCorLev is the A map for different correctness levels (of a certain report). I.e.
% if I really really want to report Attr=2, then col 2 would be
% noisyBino(corGrid(corLevN),0.2,corLevN) . If I thought that reporting Attr=2 is OKish,
% then that col in A would be noisyBino(corGrid(corLevN-1),0.2,corLevN)
% Desirability shifts will then modify these, so if des > 0 then the Attr reports which are
%              lesser than the underlying belief about Attr will be shifted towards correct, up,
%              and the >= towards wrong. Vice versa if des < 0.
if corLevN == 3
    corGrid = [0.2, 0.5, 0.8];
else
    corGrid = 1/(2*corLevN)+(0:(corLevN-1))/corLevN;
end
%  end block about wrong - OKish - correct - Neutral stuff

%% re-record for output the params ----------------------------------------------
% First, the ones that we will want fitted to data - see
%  Likelihood function e.g. spm_mdp_L_*.m , where indexP stands for
% structure of index params, as transformed to the real line.
% Also, a more extended structure with all the native-space params that
% that may be needed to construct d0, AHub, A2 etc in spm_mdp_L*.m
p4fit.attr0   = attr0;              
p4fit.dInitEv = dInitEv;        p4fit.aInitEv = aInitEv;      
p4fit.alphaPrec = alphaPrec;
p4fit.wAttr = wAttr;            p4fit.w0 = w0;    

modelStructure.indexP = nat2tr_mdp_L_xi(p4fit); % indexP in transformed space.
% So if we assume to-be-fitted to be as above,
allP = tr2nat_mdp_L_xi(modelStructure.indexP);  % i.e. back to native, for funsies, how inefficient is this way of doing it ha ha ;)
allP.lrnR = lrnR;
allP.desBias = desBias;
allP.desCorr = desCorr;      allP.desPos = desPos;
allP.AttrV0 = [0.05 0.95];
allP.resNRepo = resNRepo;    allP.resNHub = resNHub;
allP.trueAttrV = AttrVal;     % for running in generative mode etc.
% No. of levels of desirability for unpos, pos, ... indifferent/neutral :
allP.retLevN = 3;
allP.Ucor = Ucor;            allP.Upop = Upop;
allP.corLevN = corLevN;
allP.corGrid = corGrid;
allP.trialN = trialN;    allP.Tsteps1 = Tsteps1;  allP.Tsteps2 = Tsteps2;
allP.tiny = 1e-16;        % probabilities of this order are zero for all intends and purposes.
allP.noiseFloor = 0.001;  % General purpose 'small compared to other sources of noise'
allP.isOct  = is_octave();  % What, if any, Octave version we are running.

modelStructure.allP = allP;  % Store for output!

%%  ############################################################
%%  ######             Learning :   Hub MDP            #########
%%  ######     observe split & update Attr beliefs.    #########
%%  ############################################################

%%  ?? V map [ allowable policies (of depth T, so T-1 rows) ] here!
%   For the Hub, there are no actions to be taken - it's a HMM !
%   See function spm_MDP_get_T  in file spm_MDP_VB_XI.m  l. 1323
%   ----------------  so no V or U here ! ----------------------------
%   In general, policies are sequences of actions
%   (with an action for each hidden state factor, enumerated in the
%    last column. T-1 rows bec. actions not taken in last state.)
%  So of the form V(timestep, action-sequence, state-factor) I think
%  REM here stateFNHub = 2;  % No. of state factors: trueAttr, trueSI
%--------------------------------------------------------------------------

%% %% The key dynamics: B map, i.e. Transitions conditioned on action %% %
% controlled transitions: B{u}
%--------------------------------------------------------------------------
% We specify the probabilistic transitions of hidden states for each
% state (as opposed to outcome) factor. Hence actions are themselves
% factorised, state factor F being subect to actions MF :
%  REM from spm_MDP_VB_X: MDP.B{F}(NF,NF,MF)  - transitions among states
%       under MF control states. So each B{k} refers to *state factor* k,
%       which has Nk states; There must be a grand total of Mk
%       actions / control states for these, the 'pages' of B{k}, not
%       all of which might be available under V for the current time-step.
%--------------------------------------------------------------------------
% B1{1} : transitions over factor 1, trueAttr : There is only the
%        the action 'accept fate', the pt can't change this.
%        trueAttr is considered to have resNHub poss. levels:
B1{1} = eye(resNHub);

%% %%%% Learning MDP outcome probabilities: A (likelihood map) %%%%% %
%--------------------------------------------------------------------%
% REM from spm_MDP_VB_X: MDP.A{G}(O,N1,...,NF) - likelihood of O outcomes
%     given hidden states. G must stand for the *outcome factors*, whereas
%     F must stand for the *state factors*, N1 being the number of states
%     in state factor 1 etc.
% Here specify the probabilistic mapping from hidden states
% to observations (a.k.a. outcomes):
% AHub{1} : Observe positivity level reported
% AHub{2} : Actual neg or pos ret.; depends on Stage and trueAttr
% Their form is A{G}(observation, trueAttr, trueSI, level-to-report)
% See Emotion_learning_model... A map as example.
%--------------------------------------------------------------------%

% A : Actual pos ret. depends only on trueAttr and trueSI
%                 returns     trueAttr  
% AHub = zeros(retLevN-1,    resNHub    );
aHub = AHubClassLrn(allP);           % initial generative model ... from  AHubClassSerDict(allP);
aHub(1:2,:) = aHub(1:2,:) * allP.aInitEv;   %     ... likelihoods have default SIVal and AttrVal.

AHub = AHubClassLrn(allP, AttrVal);  % generative process likelihoods (from AHubClassSerDict)


%% Learning MDP priors: (utility) C1 aka goal map ---------------------------------------
%--------------------------------------------------------------------------
% Specify the prior preferences in terms of log probabilities over outcomes.
%    REM from spm_MDP_VB_X:
%    MDP.C{G}(O,T) - (log) prior preferences for outcomes (modality G), I
%    believe at each of T timesteps or moves ??
%    C factors indexed by G must corresp. to A.
%--------------------------------------------------------------------------

% Outcome is desirability of pos etc. split.
%   allP must include, for use here, just desPos
C1  = CHubLrn(allP);       % from CHubSerDict

%% Prior beliefs about initial states: d map ------------------------------
% now specify prior beliefs about initial states, in terms of counts. Here
% the hidden states are factorised into
% MAKE SURE THEY ARE COL VECTORS AS IN E.G. d{2} = [2 2]';
% state factors: stage, trueAttr, trueSI, Level2Report
%--------------------------------------------------------------------------

%  d  has 1 factor, beliefs about 'true Attr' state ofthe rater
%  allP must include, for use here: attr0, dInitEv, resNHub.
%     Optional Upop (default=1), the dispersion param. for the population distro.
d0 = dHubClassLrn(allP);                           % from dHubClassSerDict
%should produce the basic output that I should understand the physological meaning of it

%then i would need to put the trials together to have the model of the
%whole task

%%  Hub MDP Structure - to be used to generate arrays for multiple trials
%==========================================================================

%% Generative and process model structure for hub -  -  -  -  -  -
% Note conc. parameters are counts not probs. Specifying the mdp.a here
% is sufficient for the a map to undergo learning when the model is 'solved' below.

hub_a0{1} = aHub;   % hub_a0 will serve also as a working and storage copy
                    % intialise! Best to keep cell array notation even for 1 element.
hmmHub.a = hub_a0;
hmmHub.eta = lrnR;    % This is 'learning rate' (learning between trials, I hope).

% The true generative process:
hmmHub.A = AHub;
hmmHub.B = B1;                  % transition probabilities, allready cell array
hmmHub.C{1} = C1;               % preferred outcomes
hmmHub.d = d0;                  % prior over initial states
hmmHub.s = trueAttr;   % true initial state. resNHub would mean highest true Attr
hmmHub.alpha = alphaPrec;
hmmHub.T = Tsteps1 - 1;         % Why hmmHub.T = Tsteps1-1 works, rather than Tsteps1 ? hmmm ...
hmmHub.tau   = 12;
%  fix the names
hmmHub.Aname = {'returnPos'};
hmmHub.Bname = {'trueAttr'};
% The following line doesn't work properly - there may be a bug in spm_MDP_VB_trial.m
hubModalityNames =  {'ReturnValue'};
hmmHub.name.modality = hubModalityNames; % just if needed for testing - mainly see below.



%%  ############################################################
%%  ######          Actions :   REPORTING MDPs         #########
%%  ############################################################

%%  V map: allowable policies (of depth T, so T-1 rows).
%   These are sequences of actions
%   (with an action for each hidden state factor, enumerated in the
%    last column. T-1 rows bec. actions not taken in last state.)
%  So of the form V(timestep, action-sequence, state-factor) I think
%  REM here stateFNRep = 2;  % No. of state factors: beliefState, reportState
%--------------------------------------------------------------------------
V2 = nan(Tsteps2-1,resNRepo,stateFNRep);  % see above. stateFNRep relevant state factor
V2(:,:,1) = 1;  % is actions affecting state factor 1, e.g. beliefAttr, and only
                % 'acceptance' - beliefs are changed only due to Hub!
V2(:,:,2) = 1:resNRepo;  %  is actions affecting state factor 2, e.g. reportAttr

%there is reporting and learning MDP,
%what you do i neach timestep will affect what u do in the next timestep, 
%ex timestep one would be yr starting point, if u have an action MDP then u
%can take an action in timestep 1 and it will affect what happens in
%timestep 2. then the learning MDP would feel in to each timestep which
%would affect the decision.

%rule of thumb:rows are outputs and everything else is input in the maps
%ex. time, aciton etc. would be in cols not rows
%% %% The key dynamics: B map, i.e. Transitions conditioned on action %% %
% controlled transitions: B{u}
%--------------------------------------------------------------------------
% We specify the probabilistic transitions of hidden states for each
% state (as opposed to outcome) factor. Hence actions are themselves
% factorised, state factor F being subect to actions MF :
%  REM from spm_MDP_VB_X: MDP.B{F}(NF,NF,MF)  - transitions among states
%       under MF control states. So each B{k} refers to *state factor* k,
%       which has Nk states; There must be a grand total of Mk
%       actions / control states for these, the 'pages' of B{k}, not
%       all of which might be available under V for the current time-step.
%B map: takes initial state of the world and transitions it to the final
%stage of the world at that specific timestep, B map doesn't change with
%timestep. Identity maps do not change, this B map is an identity map
%--------------------------------------------------------------------------
% B2{1} : transitions over factor 1, e.g. as-believed-Attr : There is only the
%        the action 'accept fate', the pt can't change this
%        here, only the reporting of it.
B2{1} = eye(resNRepo); %ex. the transition between each timestep is the identity maps
% B{2} is the choice of reported level of Attr
B2{2} = zeros(resNRepo+1,resNRepo+1,resNRepo);
% First, from Initial state which is resolN+1(resolution of the likert scale, if it's +1 then n u did 3 then the 4th is yr initial state), valid report
% actions lead to the respective report state:
for contrS=1:resNRepo
   noisy = allP.noiseFloor * ones(resNRepo,1);
   noisy(resNRepo+1) = 0;  % Can't go to these with a valid report action
   noisy(contrS) = 1;
   noisy = noisy / sum(noisy);
   B2{2}(:,resNRepo+1,contrS) = noisy;
end
% *from* a valid report control state all actions lead to same
for contrS=1:resNRepo
   % for clarity, first from valid report states:
   B2{2}(1:resNRepo,1:resNRepo,contrS) = eye(resNRepo);
end
%run this ltr, it would add noise. ex. if u were supposed to end up in 4,
%it would intentionally give u a bit of an error maybe u could end up in 3

%% %%% outcome probabilities: A template(likelihood map) %%%%%%%%%%% %
%--------------------------------------------------------------------%
% REM from spm_MDP_VB_X: MDP.A{G}(O,N1,...,NF) - likelihood of O outcomes
%     given hidden states. G must stand for the *outcome factors*, whereas
%     F must stand for the *state factors*, N1 being the number of states
%     in state factor 1 etc.
% Here specify the probabilistic mapping from hidden states
% to observations (a.k.a. outcomes):
% A2{1} : Informational, (eg Attr)-report stage
% A2{2} : How well Intent is to be reported.
% Their form is A{G}(observation, trueState, level-to-report)
% See Emotion_learning_model... A map as example.
%curly brackets are state factors

%--------------------------------------------------------------------%

% A2{1} : Informational, (eg Attr)-report stage
%             (row)               (col)              (page)
%      observed-e.g. AttrReport  eg as-believed-Attr  eg AttrReport
A2{1} = zeros( resNRepo+1,        resNRepo,           resNRepo+1 );
% Report outcomes from initial report state, resolN+1,
%participants don't care about which timestep they are in;
%ex.if the agent believes X is going to happen they would report X, that is if
%they value being correct, they have learned from the learning MDP and then
%they report you what they have learned
% are all 'null', aka noReportMade, only coincidentally resNRepo+1:
A2{1}(resNRepo+1,:,resNRepo+1) = 1;       
% Report outcomes from valid report states are as themselves,
% whatever the true state (i.e. we observe reporting what we reported)
for toRep = 1:resNRepo
    A2{1}(toRep,:,toRep) = 1;
end
%%   A2{2} : How well Intent is to be reported.
%    Done trial-by-trial, see below.
%            (row)         (col)        (page)
%          report quality  attribBelief  attribReport
%  allP must include, for use here: desBias, resNRepo, Ucor, corLevN, resNRepo :
A2{2} = ARepClassifLrn(allP);      % from ARepClassifSerDict
% Test with e.g. squeeze(A2{2}(:,1,:)) to see correctness levels for all Hrep when Hbel is 1

%% Reporting MDP priors: (utility) C, aka goal map ------------
%--------------------------------------------------------------------------
% Finally, we have to specify the prior preferences in terms of log
% probabilities over outcomes.
%    REM from spm_MDP_VB_X:
%    MDP.C{G}(O,T) - (log) prior preferences for outcomes (modality G),
%    at each of T timesteps or moves.
%    C factors indexed by G must corresp. to A.
%
%--------------------------------------------------------------------------

% Outcome factor 1 is informational (e.g. Attr) report:
%if the states are preferred they would have a high value for being correct
%and low value for being wrong;
C2{1} = zeros(resNRepo+1,Tsteps2 );  % No preferences here!

% Outcome factor for desirability to report what one believes is best:
%   allP must include, for use here: desCorr, corLevN, Tsteps2
C2{2}  = CRepClassifLrn(allP);

%% Element, or unit, of positivity-reporting MDP :  -------

  % We will make working variables mdpAttr, which we will populated for every
  % iteration with the correct bits of hmmHub
  %% V, B, C, A, alpha, eta, omega, tau will be the same for all reporting MDPs:
  mdpPos.B = B2;       mdpPos.C = C2;
  mdpPos.V = V2;    
  mdpPos.alpha = alphaPrec;
  mdpPos.eta = allP.tiny;         mdpPos.omega = 1-allP.tiny;  % } essentialy flags.
  % Pedantic - the time constant for gradient descent in spm_MDP_VB_XI :
  mdpPos.tau = hmmHub.tau;
  mdpPos.A = A2;

  % States ------------------------------------------------------------------
  mdpPos.s = zeros(stateFNRep,1);      
  
  % Naming -------------------------------------------------------------------
  % All this needed as there's a naming bug in spm_MDP_VB_X*.m
  posModalityNames = {'whichPosiRep','positivityRepQual'};
  mdpPos.name.modality = posModalityNames;
  mdpPos.Aname = {'InfoPosRep','QualityPosRep'};
  mdpPos.Bname = {'positivityBelief','positivityReport'};

%% Store in model structure ------------------------------------------

modelStructure.hmmHub = hmmHub;
modelStructure.mdpPos = mdpPos;

modelStructure.midPPos = midPPos ;
modelStructure.d0 = d0;

%% ####################################################################
%% ###           Run the likelihood function in generative mode      ##
%% ####################################################################
Inp     = [];
Resp    = [];
details = 1;  if toPlot > 0;  details=2; end
P = spm_vec(modelStructure.indexP);    % Un-necessary to run below, but good for testing code.

% Default priors, for general use and to test likelihood function:
%                      [attr0,dInitEv,aInitEv,uPrec,wAttr,w0,lrnR, desBias]

% indexP fields:      {'attr0','dInitEv','aInitEv','alphaPrec','wAttr','w0'}  ;
modelStructure.priPar=[[1.01 ,   1.01,     1.01,          1.01,     10,  10]; ...  % A
                       [1.01 ,    2,      2,             2,        10,  10]; ...  % B
                       [ 0,      0,      0,             0,       -46  -46]; ...  % lo
                       [ 1,     100,    100,            100,      46   46]];     % hi
        %                        <- max    at 1              >   <   SD ~10  >

MDPs = spm_mdp_L_xi(P,modelStructure,Inp,Resp,details); %main structure of MDPs, likelihood function, consists of markov decision processes


%% Tidy and MDPs for output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%. Tidy up and store the Hub -------------------------------------
% There seems to be a bug naming the modalities within the above,
% fixed roughly  here:
for k=1:length(MDPs.hubHMM);  MDPs.hubHMM(k).label.modality = hubModalityNames; end
% Now store d or D in array form:
MDPs.dHub = {};                    d1 = mdp2arr(MDPs.hubHMM,'d');
MDPs.dHub{1} = [d0{1} d1{1}];      
try
   MDPs.dPos = mdp2arr(MDPs.MDPFair,'d');
catch
   MDPs.dPos = {'see_D!'};
   MDPs.DPos = mdp2arr(MDPs.MDPPos,'D');
end
% a, the likelihood concentration parameters, in array form:
MDPs.aHub = mdp2arr(MDPs.hubHMM,'a');

% Store 'Reporting MDPs' ------------------------------------------
% Separately store actions, MDP.u(F,T - 1).
%  cf. in Step_by_step* tutorial:
% DCM.Y= {MDP.u}; % include the actions made
% Rem: % MDP.u(F,T-1) - vector of actions - for each hidden factor F and for timesteps 1 ... T-1
Resp.posRep = {MDPs.MDPPos.u};
% Store observations:
% Rem: MDP.o(G,T) - matrix of outcomes - for each outcome modality G and timestep T
% G for the hub are just the fair / unfair returns.
Inp.hub = {MDPs.hubHMM.o};
% G for reports are information 'where we are' oucome and 'report quality' outcome.
Inp.posRep = {MDPs.MDPPos.o};

%% ~~~~~~~~~~~  more Plots and demos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if toPlot > 0    
   disp('Trial:');
   disp(1:trialN);
   disp('Returns delivered, 1=unfair, 2=fair:');
   R = mdp2arr(MDPs.hubHMM,'o'); disp(squeeze(R(end,end,:))');
end
% First a reminder of the learning MDP, ie. hmmHub
if toPlot > 1
%     % Core HMM ------------------------------------------------------------
%     spm_figure('GetWin','Hub (learning): First trial'); clf
%     % spm_MDP_VB_trial(MDP(1));
%     spm_MDP_VB_trial(MDPs.hubHMM(1));
% 
%     spm_figure('GetWin','Hub (learning): Last trial'); clf
%     spm_MDP_VB_trial(MDPs.hubHMM(trialN));
% end
% if toPlot > 0
    
	% fairness reporting -------------------------------------------------
    spm_figure('GetWin','positivity reporting: First trial'); clf
    % spm_MDP_VB_trial(MDP(1));
    spm_MDP_VB_trial(MDPs.MDPPos(1));

    spm_figure('GetWin','positivity reporting: Last trial'); clf
    spm_MDP_VB_trial(MDPs.MDPPos(trialN));

end


if toPlot > 2
    %% Demo of neural dynamics ------------------------------------------------
    % responses to chosen option - 1st trial, incl phase-precession

    %--------------------------------------------------------------------------
	spm_figure('GetWin','Attributing: trial 1 neural responses'); clf
	% spm_MDP_VB_LFP(mdpAttr(1),[2 3;3 3],1);
	spm_MDP_VB_LFP(MDPPos(1),[1 2;2 2],1);

	spm_figure('GetWin','Attributing: last trial neural responses'); clf
	spm_MDP_VB_LFP(MDPPos(trialN),[1 2;2 2],1);

	% illustrate phase-amplitude (theta-gamma) coupling
    %--------------------------------------------------------------------------
	spm_figure('GetWin','Attributing: all trials neural responses'); clf
	spm_MDP_VB_LFP(MDPPos(1:trialN));

end  % plotting options


return;



