function mL = spm_mdp_L_xiii(P,modStruc,Inp,Resp, repoFn,details)
% minus-log-likelihood function, developped from spm_mdp_L to go with learnReport01 - roman num index _i, _lxix etc.
% FORMAT mL = spm_mdp_L_xii(P,modStruc,Inp,Resp)
% P         - parameter VECTOR or structure, IN -INF TO INF (transformed not native) space
%             This fun. must be able to use with fminunc etc. which expect P to be vector.
% modStruc  - generative model (or sufficient ingredients thereof)
%             It is strongly advised for modStruc to include modStruc.indexP as a STRUCTURE according to 
%             which we will cast a vectorial P.
%           - Can have, ***BUT WE WILL NOT USE HERE*** modStruc.priPar [usu. has the (macro)params 
%             of the prior distros in rows of alpha, beta, lo, hi for use by dbetasc .]
% Inp aka U - inputs to gen. mod.   ('obsevations': stimuli, outcomes, where we are ...)
% Resp aka Y - observed outputs       (responses / actions performed)
% Resp, or both Resp and Inp, can be empty matrices if we want to run in generative mode. 
%
% BASIC PIPELINE:
%     . learnReport* -> formulate and practice MPD based model
%                  and assemble MAP structure / DCM structure for model-fitting to be 
%                  used within EstimParAcIn below (in this fn).
% ->  . spm_mdp_L_*     -> Likelihood fn 
%     . [experiment]Dat4AcIn*   -> bring exprerimental data to format of generative model
%     . EstimParAcIn*   -> Use fmincon (or spm_nlsi_Newton...) to fit model structure
%     . [experiment]Fit[number][letter]_[latin]_[number], e.g. myExperFit01a_xi_01 to fit 1st batch of data,
%               with param combination 'a', using likelihood function x1, gen model 01.
%     . mergeExpFits[number][letter]_[number] --> produce nice csv out of multiple fits that will have been 
%               done in parallel, including descriptive data.
%
% This dovetails with learnReport01
% This function finds the key ingredient(s) of the generative model  in modStruc
% and a given set of (transformed to -inf to inf) parameter values in P and 
% values, after adding in the observations and actions on each trial
% from (real or simulated) participant data. It then sums the
% (log-)probabilities (log-likelihood) of the participant's actions under the model when it
% includes that set of parameter values. 
%
% Fitting routines can use this function to fit model params to data. 
%%__________________________________________________________________________

% repoFn 0 is the AcIn report MDPs, 
% repoFn 1 based on full pmf sharpening (softmax-like) ... 

debugging = 0;   % will report total sum LL + logPrior if 0, or just the sum log lik if 1
try 
    details;     % details > 0 enriches mL, >= 2 plots stuff.
catch
    details=0;                     % default if actions given is to just output the log lik
    if isempty(Resp); details=1; end  % ... but actions not given, output lots of synthesized stuff.
end   

Pall   = modStruc.allP;    % comprehensive list of params to modify and 
                           % use to change d, A, C and other maps
tiny = Pall.tiny;
trialN = Pall.trialN;      resNRepo = Pall.resNRepo; 
hmmHub = modStruc.hmmHub;        
% Make sure that there is no A field if we are running in generative mode
% and inputs (pos/unfair returns etc) have been provided: 
if isfield(hmmHub,'A') && ~isempty(Inp); hmmHub = rmfield(hmmHub,'A'); end
    
mdpPos = modStruc.mdpPos;   % Reporting positivity MDPS

%% Bring parameters to native space, i.e. out of log- , logit- etc. space  
% to make suitable to pass to the model to compute the log-likelihood
%--------------------------------------------------------------------------

% REM P AND indexP have ONLY the 'index' params, e.g. to be fitted! 
%              Not the ones needed for additional specification / fixed ones! 
%              To see all the params, incld he fixed ones, list modStruc.allP .
% First, bring P from its inputted form, such as a vector, to the form that 
% the param structure modStruc.indexP also has :
if ~isstruct(P); P = spm_unvec(P,modStruc.indexP); end

% default is not to adjust the various maps unless dictated by P.
% Here record only the maps that may be meaningfully adjusted.
% Learning 'hub' :
dHubAdj=0;      aHubAdj=0;      
% CHubAdj = 0;  % Let's not adjust the learning hub C map in  spm_mdp_L_xii . 

% Reporting 'spoke(s)':
ARepAdj=0;      CRepAdj=0;      % pPosAttrAdj = 0;

% Complicated if statement to transform inputted parameters to native 
% space and replace the relevant values in Pall
field = fieldnames(P);   
Plen = length(field);   

for i = 1:Plen
    % first, log-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     if strcmp(field{i},'dInitEv')
         Pall.dInitEv = exp(P.dInitEv);       dHubAdj=1;
      elseif strcmp(field{i},'aInitEv')
         Pall.aInitEv = exp(P.aInitEv);       aHubAdj=1;
      elseif strcmp(field{i},'alphaPrec')
        Pall.alphaPrec = exp(P.alphaPrec);   
        hmmHub.alpha = Pall.alphaPrec;       
        mdpPos.alpha = Pall.alphaPrec; 
     % logit-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      elseif strcmp(field{i},'posi0')  
        Pall.posi0 = 1/(1+exp(-P.posi0));     dHubAdj=1;
        Pall.attr0 = 1 - Pall.posi0;          % the prob. of being -ve
      elseif strcmp(field{i},'mem')  %  v. simple - change directly:
        Pall.mem = 1/(1+exp(-P.mem));  
        hmmHub.omega = Pall.mem; 
    % and scaled-logit transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     elseif strcmp(field{i},'desBias')
%         error('P should not contain desBias field at this stage');
%         % Pall.desBias= -1+ 2/(1+exp(-P.desBias));   ARepAdj=1;  
    % Untransformed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      elseif strcmp(field{i},'wAttr')
        % Here we check before changing, bec. if we don't 
        % need to change we can save an expensive function call later.
        if abs(Pall.wAttr - P.wAttr) >  tiny
           Pall.wAttr = P.wAttr;   aHubAdj=1;     pPosAttrAdj = 1;
        end
      elseif strcmp(field{i},'wp0')  
        if abs(Pall.wp0 - P.wp0) >  tiny
         Pall.wp0 = P.wp0;    Pall.w0 = -P.wp0;    aHubAdj=1;    pPosAttrAdj = 1;
        end
      else
       error([field{i} ' not catered for :(']);
    end
end

%% ** Assemble and run the generative model with correct stimuli and returns  **

%% Extract the single-trial MDPs from modStruc, components of the overall model
%  that we will need to put together. Also extract other auxiliaries.

% Adjust component MDP maps acc to parameters as necessary  -------------
% NOT NEEDED - useful only if multiple attributes
% if pPosAttrAdj 
%     % pPosHiLoAttr provides all the pos return probabilities for all the values of Attr
%     % states that the participant considers if Attr can be either high or low, and 
%     % these low and high- attribute states can only take values from attrV below. iLo and iHi
%     % are the indices over attrV for each row of pPosAttr, for Attr responsible acc. to attrV
%     % % legacy from one used in allP in [mdpHI, mdpSI] = hub2spoke_dAttr(mdpHI,mdpSI, mdpHubElement, allP)
%     Pall.pPosHiLoAttr = pPosHiLoAttr([Pall.w0,Pall.wAttr]);  % provides pPosAttr, attrV, iLo, iHi, iAttr
% end

if   aHubAdj
    a0 =  AHubClassLrnPos(Pall);  
    a0(1:2,:) = a0(1:2,:) * Pall.aInitEv; 
    hmmHub.a{1} = a0; 
    if isempty(Inp)     % i.e. we have not been given returns etc, so need to generate them
        error('Inp, i.e. the ratings, should be provided to the SEL likelihood function');
    end
end

% Let's not adjust the learning hub C map in version spm_mdp_L_xii ... 
% if   CHubAdj; hmmHub.C{1}  = CHubLrn(Pall); end 

if  dHubAdj 
    d0 = dHubClassLrnPos(Pall);   % Will need copy before it's updated!
    hmmHub.d     =  d0; 
end

% Adjust report MDP maps acc to parameters, if necessary  --
if  ARepAdj;      mdpPos.A{2} = ARepClassifLrn(Pall);    end
if CRepAdj;       mdpPos.C{2} =  CRepClassifLrn(Pall);   end


%%  1. 'Hub', or Learning MDP ---------------------------------------------------------

% Assemble trials and add observations and actions ------------------------------------
[hmmHub(1:trialN)] = deal(hmmHub);          % Create MDP with specified number of trials

% Cater for deliberately missing function arguments, when fn. used in 
% generative mode : 
if ~isempty(Inp)
   % Note that in SEL, the observations of the learning hub
   % DO NOT DEPEND ON THE DECISIONS OF THE AGENT, so they should be
   % copied straight from the input. 
   [hmmHub.o]      = deal(Inp.hub{1:trialN});     % Add observations in each trial
end

%--------------------------------------------------------------------------
%  Run model for hub, so as to map beliefs about each available state
%  in each trial, given parameters and stimuli :  
hmmHub   = spm_MDP_VB_XI(hmmHub); % run model with given parameter values
%--------------------------------------------------------------------------

% Store true generated states and flag appropriately if in generative mode:
if isempty(Resp)
    Resp.sHub = {hmmHub.s};     
    genMode = 1;
else
    genMode = 0;
end
    

%%  Reporting 'spoke(s)' --------------------------------------------------------------

LPos=0; % initialise summed (log) probability of reporting MDP actions

for tr = 1:trialN
    
  
  % Form the D map, i.e. belief distibution over of the positivity reporting. We must use
  % the posterior of the PREVIOUS trial, or the very starting one:
  if tr > 1
     % mdpF.d = hub2spoke_dPos(hmmHub(tr-1).a{1}, hmmHub(tr-1).d, Pall);
     Dpos = hub2spoke_D_Pos(hmmHub(tr-1).a{1}, hmmHub(tr-1).d, Pall);
  else
     Dpos = hub2spoke_D_Pos(a0, d0, Pall);  % use copies from before hmmHub model ran at all!
  end
  
  %%  Sharpened probabilities / Softmax like Reporting 'spoke(s)' --------------------------------
  if repoFn == 1
      polP  = Dpos{1} .^ Pall.alphaPrec;   % as if lnP is value ;)
      polP  = polP / sum(polP);
      
      if genMode % if in generative mode, sample the Report:
          u = pBinSample(polP); 
          Resp.posRep{tr} = [1 u]';   % The '[1 ' is for compatibility with repoFn = 0 simulations.
      else  % if not in generative mode, the Report must be given:
          u = Resp.posRep{tr}(2,1); 
      end
      l1 = log( polP(u) + 1e-16);
 
      LPos = LPos + l1; 

      if details  % Store progress in detail for debugging etc. ~~~~~~~~~~
            mL.lPos(tr) = l1;                   
      end % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  elseif repoFn == 0
      %%  AcIn MDP Reporting 'spoke(s)' --------------------------------------------------------------
      % Copy from updated-parameter element prototypes:
      mdpP = mdpPos; 

      % Set the crucial prior beliefs that we will try to report, w. noise, bias etc.:
      mdpP.D = Dpos; 

      % Observations / stimuli in this model include observations of one's
      % own responses, as well as pos or unfair returns. Only the third row
      % of observations of the 'learning MDP' contains the pos/unfair ret.,
      % and this does not depend on pts. beliefs and actions, hopefully!
      % Should always go neutral - pos_or_unfair 

      % Accumulate log-likelihood for Reporting MDPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      % REM: MDP.P(M1,...,MF,T)    - probability of emitting action M1,.. over time

      % Set true states for the attribution MDPs. REM columns are within-trial timepoints, 
      % rows are state factors. So top row is AttrTrue, etc.
      mdpP.s(1,1)   = hmmHub(tr).s(1,1);       
      mdpP.s(2,1)   = resNRepo+1;              

      % Set the crucial actions: 
      if ~genMode
          mdpP.u = Resp.posRep{tr};
      end

      %% ~~~~~~~~~~~~~~~ Solving and Storing the Attibuting MDPs ~~~~~~~~~~~~~~~~~~~~~

      % Now solve and store:
       mdpP  = spm_MDP_VB_XI(mdpP); 

      if genMode || details % i.e., if we are in generative mode or user wants gory details
        if genMode % test again, lol
          Resp.posRep{tr}  = mdpP.u;      %#ok<*SAGROW>
        end
        if tr == 1
          MDPPos(1:trialN) = deal(mdpP);  % }  store first trial and do memory prealloc.
        else
          MDPPos(tr) = mdpP; %#ok<*SAGROW>
        end
      end

      % clear mdpP; % inefficient but safe ...

      % Get probability of  action to be taken for reporting MDP:
      l1 = log( mdpP.P(1, Resp.posRep{tr}(2,1)) + 1e-16);  

      LPos = LPos + l1; 

      if details  % Store progress in detail for debugging etc. ~~~~~~~~~~
            mL.lPos(tr) = l1;                   
      end % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  end  % if ... elseif ... selectiof of different report functions. 
end           % end loop over trials

if details == 0                   % just the sum log lik ...
    mL = -LPos ;
elseif details <=2  % Detailed output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % record stuff
    mL.sum = LPos ;                       %  plain sum log likelihood   
    mL.hubHMM=hmmHub;   
    if repoFn == 0; mL.MDPPos= MDPPos; else  mL.MDPPos = []; end   % Not if it hasn't been used ...
    mL.P = P;                  mL.Pall = Pall; 
    mL.Inp.hub = {hmmHub.o};   mL.Resp = Resp;   % Copies handy if we ran in generative mode.
    mL.a0Hub = a0; 
    mL.d0Hub = d0;
else
    error(['Details option ' num2str(details) ' appears not catered for']);
end
    % plot stuff if asked:
if details >= 2 && repoFn == 0
      
	    spm_figure('GetWin','Positivity reporting: All, trial-by-trial'); 
	    spm_MDP_VB_game(MDPPos);
    
        % aHubProbDisp : display how the likelihoods in a are learnt over trials.
        aHub = mdp2arr(hmmHub,'a');  aHub = aHub{1}; 
        aHubProbDisp2
         
        disp(' ');
        disp('~~~~~~~~~ Key to muliple trial plots: ~~~~~~~~~~~')
        disp('Topmost  : Circles: true init. state along each factor, from mdp.s. Colours are codes:')
        disp('  Modulo numeric <-> {''r.'',''g.'',''b.'',''c.'',''m.'',''k.''}; ')
        disp('  Rows of circles corresp. to factors, first factor TOP-MOST, ')
        disp('  e.g. if initial positivityLevelReportState=5 -> magenta')
        disp('2nd down : Observations / outcomes for each factor, from mdp.o ')
        disp('  TOP-MOST is LAST factor here :(, positivity of outcome, 3=null/not present in timestep, 1=NEGATIVE')
         % disp('           ')
        disp('3rd down : ERP ;     4th down : DA ')
        disp('5th down : Post. beliefs about states, (D / updated d). y axes bizarelly seem to go HIGHER DOWN! ')
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')  


end
    

% A bit of tidying up
if repoFn ==0; clear('hubHMM','MDPPos');  else clear('hubHMM'); end

return;

