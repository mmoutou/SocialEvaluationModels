function optPetc = SEL2_Grid04a_xiii(selAIFile, batch2fit, initTrParFile, grN, grSc, batchSize, fitDir, codeTest )
%  optPetc = SEL2_Grid04a_xiii(datFile, batch2fit, initTrParFile, grN, grSc, batchSize, fitDir, codeTest ))   :
%                              MAP fitting, with weakly informative priors 
%  Fitting of YSELT data, Carlisi et al, over a grid, supplemented with key descriptives.
%  
% SEL2_Grid04*_**  MAP fitting, with weakly informative priors 
%                  initTrParFile has initial vals of params in fitting (transf. space, as per
%                  [numID,'posiSelf','posiOther','dEvSelf','dEvOther','aEvSelf','aEvOther','alphaPrec','genLR','repLR','wp0','wAttr','mem'... ]  
%                         BIC derived fairly rigorously from iBIC approx:
%                 BIC_kc := -2 ln Lopt + (k-c) ln(n/2*pi) + (c/N) ln (nN)
%                         = -2 ln Lopt + (k-c)*(ln n - ln(2*pi)) + (c/N) ln (nN)
%  AIC_kc derived rather more roughly:
%                 AIC_kc := -2 ln Lopt + 2*((k-c)  + c/N)
%                         k being the number of params fitted for each pt, and c the number of
%                         params derived (supposedly optimally ...) for the whole sample. 
%                        Fitting being done with learning from one block to the other,
%                        adding admixture of the posterior to the prior.
%                        Estimate over a grid. xiii is the version of the likelihood fn used.
% To add / remove params to grid over, change: 0. Word search 'Fixed value' and set *TrRef values
% 1. parNOTTOFIT 2. par2fitN  3. parOrdGrid 4. fitFileName 5. If need be, if par2grid <= par2fitN-1... statementS
% NB hd and totParN NOT to be changed for this purpose.
% 
% Example use with synthetic data, 
% from .../SocEvalModels_(not_git)_rough_work/dataFits/fittingTests/synRepDictCIs
%   fitD = cd; % '/home/michael/googledirs/MM/SocEvalModels_(not_git)_rough_work/dataFits/fittingTests/synRepDictCIs/'
%   synParF = 'ldat_gen_100_synRepDictCIs.mat';
%   optPtest = SEL2_Grid04a_xiii('dat100_synRepDictCIs.mat',1,synParF,12,5,20,fitD,0)

% BASIC PIPELINE:
%     . learnReport* -> formulate and practice MPD based model
%                  and assemble MAP structure / DCM structure for model-fitting to be 
%                  (may be used within EstimParAcIn below (in this fn)).
%     . spm_mdp_L_*     -> Likelihood fn 
%     . [experiment]Dat4AcIn*   -> bring exprerimental data to format of generative model
%            e.g. ...\Dropbox\BASOR\AcInSOR\AcInRepeatedDictator\attrssrib\attrSSRIbDat4AcIn09.m
%-->  [ . [expt]Grid[number][let]_[roman]_[pts batch2fit}  -> Grid fit INCL learning from one block to the next
%  or { . EstimParAcIn*   -> Use fmincon (or spm_nlsi_Newton...) to fit model structure
%     { . [experiment]Fit[number][letter]_[latin]_[number], e.g. attrssriaFit07a_iv_01 to fit 1st batch of data,
%               with param combination 'a', using likelihood function iv, gen model 07.
%     . mergeExpFits[number][letter]_[number] --> produce nice csv out of multiple fits that will have been 
%               done in parallel, including descriptive data.
fitStr = '04a';   % short string to identify the fitting version in outputs / file


%% Cater for Fixed value parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  A parameter vector in native space:
%          1           2          3        4           5        6         7            8      9      10     11     12    
parHd=  {'posiSelf','posiOther','dEvSelf','dEvOther','aEvSelf','aEvOther','alphaPrec','genLR','repLR','wp0','wAttr','mem'};
%          13   14   15    16 
%eg par = [0.75,    0.5,         1,        1,         2,       2,        5,          0.1,    0.5,    0,     6,   0.9999 ]; 
fitHd = { 'LL','LP','AIC','BIC'};         hd = {parHd{1:end}, fitHd{1:end}};
hd4outp = hd;   % we'll need a backup ;)

parNOTTOFIT = [4, 6,    5, 8, 12];  % 5,8, even 12 will have to be freed for nice fits :)  
groupParN = 0;  totPtN = 100;  % Number of parameters, a subset of parNOTTOFIT, which have been
    % from a whole group of totPtN worth data (e.g. a median of a detailed fit) and hence fixed here.
    % Needed for BIC and AIC corrections.
if groupParN == 0
   warning('Fit set for NO group-level fixed perameters');
else
   warning(['groupParN set to ' num2str(groupParN), ' based on ' num2str(totPtN) ' (notional even) participants']);
end

% Fixed, rather than yoked, params (4 and 6 will be yoked to 3 and 5 for initial fits) :
aEvSelfTrRef = log(2);               % 5.  Fairly flexible default.
genLRTrRef = logit(0.1);             % 8.  default lowish learning from block to block
memTrRef = 6.91;                     % 12. default negligible forgetting from trial to trial.
yokePar  = [[3 4];  [5 6]];          % 3. determines 4.,  5. determines 6. 

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
try
 load(initTrParFile ); % to povide ldat, which has initial values of the parameter of the
                      % fit, possibly derived from previous attempts at fitting.
                      % It may over-write hd, so fix this below.
 hdLdat = hd;         hd = hd4outp; 
catch
 ldat = [];
end

try codeTest; catch codeTest=0; end  % affects grid size, batch2fit size etc.
try grN;  catch, grN=12;    end
try grSc; catch, grSc = 5;  end      % will try grids of: par0 + (-grN:grN)*(abs(par0)+1)/grSc ; 
if codeTest
   warning('running with code testing settings; For debug, try grN=1 and grSc either 2 or 200')
end

% NB each pt has 4 partners !! So batchSize = 3 means 12 fits,
try batchSize;  catch,   batchSize = 18;  end
if codeTest; batchSize = 2; end

%% Set working directory and load data to fit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cwd = pwd;
try fitDir; catch, fitDir = []; end
if isempty(fitDir)
    if ~codeTest
        fitDir = 'C:\Users\mmpsy\Nextcloud\MM\googledirs\SocEvalModels_(not_git)_rough_work\dataFits';
    else
        fitDir = 'C:\Users\mmpsy\Nextcloud\MM\googledirs\SocEvalModels_(not_git)_rough_work\dataFits\fittingTests';
    end
    warning(['working dir set to ' fitDir]);
end
cd(fitDir); 

% selAI file provides selD which has the YSELT data in format for active inference modelling          
load(selAIFile);        % has selD. Data must be (re-formated) to fit w model Inp, Resp etc.
modStrucFile = 'modStruc8a.mat';           load(modStrucFile);     % This should contain modStruc8
% expDatFile = ['selD03_dat2.mat'];               
% load(expDatFile);       % additional parts of the data, which don't depend on the generative model.
ptNloaded = length(selD); 
totBlN = length(selD{1}.Inp);    
trNPerBlock =  modStruc8.allP.trialN;
datNPerPt = totBlN * trNPerBlock;  % total data points per participant

startPt = (batch2fit-1)*batchSize+1; 
if startPt > ptNloaded; error('Starting pt to fit is beyond length of data'); end
endPt =  batch2fit*batchSize;
if endPt > ptNloaded
    endPt = ptNloaded; 
    warning('End pt to fit was past the last one in the data, so set back to the last one.'); 
end
todo = startPt:endPt; 

% Name of the output file ---------------------------------------------------------------
fitFileName = ['sel_Grid' fitStr '_' num2str(todo(1)) 'to' num2str(todo(end))];
if grN ~= 12 || grSc ~= 5    % if not the defaults
    fitFileName = [fitFileName '_grSize' num2str(2*grN+1) 'grSc' num2str(grSc)];
end

expGrid = {};
indexPLen = length(fields(modStruc8.indexP)); 
par2fitN=length(parHd) - length(parNOTTOFIT);  % num of free params to fit at pt level. 
totParN=length(parHd); % total number of params in indexP... to check ... 

% 'reminder' of 'Likert scale' type bins for reported beliefs of the expected posiSelf and posiOther:
resNRepo = modStruc8.allP.resNRepo; 
% midGrid = 1/(2*resNRepo)+(0:(resNRepo-1))/resNRepo;
% To help store descriptives, the order in which conditions will be stored:
ord2store = [111,1; 132,2; 110,3; 120,4; 130,5; 211,6; 232,7; 210,8; 220,9; 230,10];
measN = 3* size(ord2store,1);  % one for 'early', 'mid', 'late' part of each block.
hdLong = {'numID'};
hdLong((end+1):(end+length(hd4outp))) = hd4outp; 
hdLong((end+1):(end+measN))={'eaS1neg','eaS1pos','eaSneg','eaSneu','eaSpos',...
                             'eaO1neg','eaO1pos','eaOneg','eaOneu','eaOpos',...
                          'midS1neg','midS1pos','midSneg','midSneu','midSpos',...
                          'midO1neg','midO1pos','midOneg','midOneu','midOpos',...
                       'ltS1neg','ltS1pos','ltSneg','ltSneu','ltSpos',...
                       'ltO1neg','ltO1pos','ltOneg','ltOneu','ltOpos' };
optPetc = nan(length(todo), totParN+4+measN);   % for optimised parameters AND 4 measures of model fit,
         % LL, MAPd, AIC, BIC incl the non-MDP params and the ones we won't explore here.
optNatP = optPetc;   % copy with params in native space
numID = nan(length(todo),1);                   % just to hold  numIDs

%% Main loops over participants ---------------------------------------------------------

%% Descriptives -------------------------------------------------------------------------
ptCnt = 0;        % counts participans
for ptN = todo
    ptCnt = ptCnt+1; 
    try % try setting a numerical ID
        numID(ptCnt) = selD{ptN}.numID;    
    catch
        numID(ptCnt) = ptN;
    end

     %% Clear the decks and set stuff common to ALL blocks 
     mapMod.mdpStruc = modStruc8;
     mapMod.mdpStruc.indexP = nat2tr_SEL_i([]); 
     mapMod.Inp = selD{ptN}.Inp;
     mapMod.Resp = selD{ptN}.Resp;

     % Derive the descriptives. For each of 
     % self-veRepeated, self+Repeated, self-_unrep, self_neut_unrep self+_unrep , other similar ,
     % derive early, mid and late averages using weights linearly decreasing away from
     % the defining point of each section. Conditions are, in order and if present,
     % 111,132,110,120,130, 211,232,210,220,230
     % and the corresp. data are the mapMod.Resp{block}.posRep{trial}(2) 
     for blkN = 1:totBlN
         blTrN = length(selD{ptN}.Inp{blkN}.posRep);  
         d = nan(1, blTrN);
         for trN = 1:blTrN; d(trN) = mapMod.Resp{blkN}.posRep{trN}(2); end

         ea = tent1dFilter(d,1,round(blTrN/2));
         mid = tent1dFilter(d,round(blTrN/2),round(blTrN/2));
         lt = tent1dFilter(d,length(d),round(blTrN/2));

         % find where to store these:
         cond = mapMod.Inp{blkN}.condCode;
         col  = find(ord2store(:,1)==cond);
         optPetc(ptCnt,totParN+4 + col) = ea;
         optPetc(ptCnt,totParN+4 + col + size(ord2store,1)) = mid;
         optPetc(ptCnt,totParN+4 + col + 2*size(ord2store,1)) = lt;         
     end

     % Store Key Descriptives and IDs, with ugly backup for legacy compatibility:
     ldatBak = ldat;           hdBak=hd; 
     ldat = [numID optPetc];   hd = hdLong;   % THEN RESORE ... urgh ugly legacy stuff ...
     save([fitFileName '_ldat.mat'],'ldat','hd'); 
     mat2csv2Dfl(ldat,[fitFileName '_ldat.csv'],0,1,hdLong);
     ldat = ldatBak;           hd = hdBak;  % urgh ugly ...
     
end  % loop over descriptives 

%% Grid fitting -------------------------------------------------------------------------
ptCnt = 0;        % counts participans
for ptN = todo
    ptCnt = ptCnt+1; 
    try % try setting a numerical ID
        numID(ptCnt) = selD{ptN}.numID;    
    catch
        numID(ptCnt) = ptN;
    end
    % REM    parNOTTOFIT = [4, 6,    5, 8, 12];  % I guess that 3,8,12 will have to be freed for nice fits :)  
    % REM        1           2          3        4           5        6         7            8      9      10     11     12    
    % parHd={'posiSelf','posiOther','dEvSelf','dEvOther','aEvSelf','aEvOther','alphaPrec','genLR','repLR','wp0','wAttr','mem'};
     parOrdGrid = [1 2 3 7 9 10 11 1 2 3 7 9 10 11 ];  % in order above in hd. Cycle twice!
    
     if codeTest; parOrdGrid = [1 3 10 1]; end    

     %% Clear the decks and set stuff common to ALL blocks 
     mapMod.mdpStruc = modStruc8;
     mapMod.mdpStruc.indexP = nat2tr_SEL_i([]); 
     mapMod.Inp = selD{ptN}.Inp;
     mapMod.Resp = selD{ptN}.Resp;

% %      Derive some descriptives. For each of 
% %      self-veRepeated, self+Repeated, self-_unrep, self_neut_unrep self+_unrep , other similar ,
% %      derived linearly weighed early and late averages. 
% %      these are, in order and if present, 111,132,110,120,130, 211,232,210,220,230
% %      and the corresp. data are the mapMod.Resp{block}.posRep{trial}(2) 
%      for blkN = 1:totBlN
%          blTrN = length(selD{ptN}.Inp{blkN}.posRep);  
%          d = nan(1, blTrN);
%          for trN = 1:blTrN; d(trN) = mapMod.Resp{blkN}.posRep{trN}(2); end
% 
%          Commented block before using dedicated filter function for local averages:
%          Weights to derive early and late averages of trial:
%          z = sum(round(blTrN/2 + 0.1):-1:1);  % normalizing constant needed below
%          eaW = [round(blTrN/2 + 0.1):-1:1 zeros(1, 20 - round(blTrN/2 + 0.1))];
%          eaW = eaW / z;
%          ltW  = fliplr(eaW);
%           ea = sum( d .* eaW);
%          lt = sum( d .* ltW); 
% 
%          ea = tent1dFilter(d,1,round(blTrN/2));
%          mid = tent1dFilter(d,round(blTrN/2),round(blTrN/2));
%          lt = tent1dFilter(d,1,round(blTrN/2));
% 
%          find where to store these:
%          cond = mapMod.Inp{blkN}.condCode;
%          col  = find(ord2store(:,1)==cond);
%          optPetc(ptCnt,totParN+4 + col) = ea;
%          optPetc(ptCnt,totParN+4 + col + size(ord2store,1)) = mid;
%          optPetc(ptCnt,totParN+4 + col + 2*size(ord2store,1)) = lt;         
%      end

     
     %% Serial, 'factorized' grid fitting over parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     parCnt = 0; 
     %% Loop over one parameter at a time ````````````````````````````````````````````````
     for par2grid = parOrdGrid
        parCnt = parCnt + 1;
        expGrid{ptN,par2grid}.parPDens = nan(3,2*grN+1);  % to work on finding MAP etc. Rows of 
                                                          % param value, fit measure, sum LL.
   
        %% First time aroune, initialize param values for the (iterative convergence) fit. 
        %  Some param. values are just derived from 1st block.
        %  Use the average values of potitivity as starting points for the priors. NOTE
        %  THESE ARE NOT PRIORS IN THE COG MODEL, ONLY INITIAL VALUES FOR FITTING.
        %  Also master copy of MAP model specification.
        if parCnt == 1
           natPinit =  nat2tr_SEL_i([]); %   tr2nat_mdp_L_xii(modStruc8.indexP);  % intialize native space params
           try 
               mPosi = selD{ptN}.Resp{1}.posRep{1}(2); 
               for tr=2:trNPerBlock; mPosi = mPosi + selD{ptN}.Resp{1}.posRep{tr}(2);  end
               mPosi = mPosi/trNPerBlock;
               natPinit.posiSelf  = 1/12 + (mPosi-1)/6;    % rescaled from 1-6 to 0-1, very roughly ... 
               natPinit.posiOther = 1/12 + (mPosi-1)/6;
           catch  
               warning('Could not set initial values from data - BEWARE!');
               natPinit.posiSelf = 0.667;    natPinit.posiOther = 0.667;  
           end

           % Now fill in the rest of the initial conditions for the fit, or 
           % if initial cond. file has been given well, use the values from that.
           if isempty(ldat)    % if no good init. cond. loaded
                % params 3 to 10 - used inside MDPs :
                natPinit.dEvSelf = 1;      % 
                natPinit.dEvOther = 1;     % 
                natPinit.aEvSelf  =2;      % 
                natPinit.aEvOther =2;      % 
                natPinit.alphaPrec =1;     % 1: fairly weak value ...
                natPinit.genLR  = 0.1;     % 
                natPinit.repLR  = 0.5;     % 
                natPinit.wp0  = 0.2;      % No Eye Deer if any good ...
                natPinit.wAttr  = 6;       % strong value ..
                natPinit.mem  =0.9999;     % no forgetting default
                 
                % mapMod.indexP is used directly by spm_mdp_L_vi to spm_unvec its 
                % parameter vector argument, whereas fitting routine also uses mapMod.iniTrP, the 
                % initial values for the fit in transformed space :
                mapMod.mdpStruc.indexP = nat2tr_SEL_i(natPinit); 

           else    % if good init cond. loaded
                if ~strcmp(hdLdat{9},'genLR'); error('hd{9} ~= genLR'); end
                % mapMod.mdpStruc.indexP = modStruc8.indexP;  % THIS IS A FILLER CONSTRUCTOR
                
                for k=1:totParN
                    mapMod.mdpStruc.indexP.(hdLdat{1+k}) = ldat(ptN,1+k);
                end

                % fixed values derived from whole group prev. fits or other params
                %                1           2          3        4           5        6         7            8      9      10     11     12    
                %   parHd=  {'posiSelf','posiOther','dEvSelf','dEvOther','aEvSelf','aEvOther','alphaPrec','genLR','repLR','wp0','wAttr','mem'};
                % 
                % parNOTTOFIT = [4, 6,    5, 9, 12];  % I guess that 3,9,12 will have to be freed for nice fits :)  
                % % Fixed, rather than yoked, params (4 and 6 will be yoked to 3 and 5 for initial fits)
                % aEvSelfTrRef = log(2);               % 5.  Fairly flexible default.
                % genLRTrRef = invlogit(0.1);          % 8.  default lowish learning from block to block
                % memTrRef = 8.5;                      % 12. default negligible forgetting from trial to trial.
                mapMod.mdpStruc.indexP.dEvOther =  mapMod.mdpStruc.indexP.dEvSelf;
                mapMod.mdpStruc.indexP.aEvSelf  =  aEvSelfTrRef;
                mapMod.mdpStruc.indexP.aEvOther =  mapMod.mdpStruc.indexP.aEvSelf;
                mapMod.mdpStruc.indexP.genLR    =  genLRTrRef;
                mapMod.mdpStruc.indexP.mem      =  memTrRef;   
           end

           fields4MDPs = fieldnames(mapMod.mdpStruc.indexP); 
           refMapMod = mapMod;     % master copy
    
        end   % if parCnt == 1 , i.e. to initialize parameter values for the grid fitting.    
        
        %% Key loop: Iterate over the 1-D grid of the current parameter ''''''''''''''''''''''' 
        for iGr = -grN:grN      
            mapMod = refMapMod; % REM mapMod is max a posteriori focused model whose elements
                                % will be modified acc. to learning, coalitional bias etc. and
                                % then be directly passed to likelihood (incl. MAP density) function.

            %% Parameter value acc. to grid, but before adjustment for which block we're in
            if ~isempty(find(par2grid == parNOTTOFIT, 1))  % Shouldn't get to this!
                disp('  OUCH '); 
                error(['Why did we get to par2grid=' num2str(par2grid) ' ?']);
            end
            refp = refMapMod.mdpStruc.indexP.(fields4MDPs{par2grid}); 
            thisp =  refp + iGr * (abs(refp) + 1)/grSc ; 
            mapMod.mdpStruc.indexP.(fields4MDPs{par2grid}) = thisp; 
            parName = [fields4MDPs{par2grid} 'Tr'];

            %%  ################ Fixed values ###############
            %   See above for full comments.
            mapMod.mdpStruc.indexP.dEvOther =  mapMod.mdpStruc.indexP.dEvSelf;
            mapMod.mdpStruc.indexP.aEvSelf  =  aEvSelfTrRef;
            mapMod.mdpStruc.indexP.aEvOther =  mapMod.mdpStruc.indexP.aEvSelf;
            mapMod.mdpStruc.indexP.genLR    =  genLRTrRef;
            mapMod.mdpStruc.indexP.mem      =  memTrRef;   
            
            %%
            if iGr == -grN && parCnt ==1  % initialise optimal param estimate - elements will be replaced 
                % as we optimise. Last two are for the fit measures. 
                % REM fitHd = {                                  'LL','LP','AIC','BIC'}; 
               optPetc(ptCnt,(1:totParN+4))= [spm_vec(mapMod.mdpStruc.indexP)' nan(1,4)];
            end
                
            % display progress:
            t = clock; tstr = num2str(t(2:5)); 
            if ptCnt == 1; tStart = tstr; end
            disp(['Now pt.=' num2str(ptN) ' for ' fitFileName ' at ' tstr]);

            % Extra data isn't quite ready, so try for future reference ... 
            try 
               mapMod.file = expD{ptN}.file;
               mapMod.identiCode = selD{ptN}.indentiCode;
            catch
               % just carry on
            end
               
            %% ~~~~~~~~~~~ The crucial Log Lik Density ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            %   ================================================================================
            pTr = spm_vec(mapMod.mdpStruc.indexP)';
            mL = sel8bl03(pTr,mapMod.Inp,mapMod.Resp,mapMod.mdpStruc,...
                         mapMod.mdpStruc.repoFun,1);  % last entry is 1 to give a bit of detail.
            expGrid{ptN,par2grid}.parPDens(1,iGr+grN+1) = thisp;                         
            expGrid{ptN,par2grid}.parPDens(2,iGr+grN+1) = mL.pDens;
            expGrid{ptN,par2grid}.parPDens(3,iGr+grN+1) = mL.sLL;
            %   ================================================================================
            
           % Show gory details if testing code, otherwise beginning, middle and end of this grid
           if  iGr == -grN || iGr== 1 || iGr == grN || codeTest
               disp(['Tranf.param. val.: ' num2str(thisp) '; sum post. dens.: ' num2str(mL.pDens) '; sum LL: ' num2str(mL.sLL)]);
           end
                         
        end  % end for iGr = -grN:grN
          
        %% ~~~~~~~~~~~~~ Find optima and Save as we go along ~~~~~~~~~~~~~~~~ 
        [MAPv, MAPi] = max(expGrid{ptN,par2grid}.parPDens(2,:)); 
        optp = expGrid{ptN,par2grid}.parPDens(1,MAPi);
        if par2grid <= totParN
           optPetc(ptCnt, par2grid)  = optp;
           refMapMod.mdpStruc.indexP.(fields4MDPs{par2grid}) = optp; 
        else  % this should be NOT pocBTr
           error('Why did we get here??');      
        end
        % fixed params: also record the yoked ones ---------------------
        yok = find(yokePar(:,1)==par2grid) ;  % see if the iterated param has others yoked to it
        if ~isempty(yok)
            for k=1:length(yok)
                yoked = yokePar(yok(k),2);
                optPetc(ptCnt,yoked) = optp;
                refMapMod.mdpStruc.indexP.(fields4MDPs{yoked}) = optp; 
            end
        end
        
        % New best param fit measures ----------------------------------
        optPetc(ptCnt,totParN+2) = MAPv;   % the fit measures, as per hd = {... 'LL','LP','AIC','BIC'}; First is MAP density, then ...
        sumLL = expGrid{ptN,par2grid}.parPDens(3,MAPi); % sum LL (make not to calc BIC etc below), then,
        optPetc(ptCnt,totParN+1) = sumLL;
        %  AIC_kc derived a little slopily :) c=1 as only wS is a fixed value param, to group median:
        %         AIC_kc  := -2 ln Lopt + 2*((k-c)  + c/N)
        optPetc(ptCnt,totParN+3) = -2*sumLL + 2*par2fitN + ...
                                             2*groupParN/totPtN ;     % roughly corrected Akaike IC (AIC)
        % BIC derived fairly rigorously from iBIC approx:
        %                 BIC_kc := -2 ln Lopt + (k-c) ln(n/2*pi) + (c/N) ln (nN)
        %                         = -2 ln Lopt + (k-c)*(ln n - ln(2*pi)) + (c/N) ln (nN)
        optPetc(ptCnt,totParN+4) = -2*sumLL + par2fitN*(log(datNPerPt) - log(2*pi)) + ...
                                   + (groupParN/totPtN)* log( totPtN * datNPerPt) ;  % Bayesian IC
                                   % incl. terms for BIC for small data sample per pt, for good measure
                                   % as some params only estimated from certain data blocks (prob. unnecessary ...)
        % New best params in native space:
        optNatP(ptCnt,totParN+1:totParN+4) = optPetc(ptCnt,totParN+1:totParN+4); 
        optNatP(ptCnt,1:totParN) = tr2nat_SEL_i(optPetc(ptCnt,1:totParN)); 
        nDone = (ptCnt - 1)*length(parOrdGrid) + parCnt;
        nTot  = length(parOrdGrid)*length(todo);
        disp(['Done ' num2str(nDone) ' out of ' num2str(nTot) ' since ' tStart ]);                                                             
        disp(['Best tr. ' parName ' so far: ' num2str(optp) ' MAP density:' num2str(MAPv)]);
        disp([' Opt. native par so far: ' num2str(optNatP(ptCnt,:))]);
        
        save([fitFileName '.mat'],'expGrid','optPetc','optNatP','hd','hdLong','numID');
    end

    % Version with params in transformed space, also the basis for further fits:
    optPar = [numID optPetc(:,1:totParN+4)];
    hdBak=hd; hd = {'ptID', hd4outp{1:end}};  % ugly backup for legacy compatibility
    mat2csv2Dfl(optPar,[fitFileName '.csv'],0,1,hd);
    % Version with params in native space, for intuitive perusal:
    optPar = [numID optNatP(:,1:totParN+4)];
    mat2csv2Dfl(optPar,[fitFileName '_natP.csv'],0,1,hd);
    hd = hdBak;   % IMPORTANT TO RESTORE.

    % Now also store versions with params in transformed space, but also 
    % Key Descriptives and IDs, with ugly backup for legacy compatibility:
    ldatBak = ldat;           hdBak=hd; 
    ldat = [numID optPetc];   hd = hdLong;   % THEN RESORE ... urgh ugly legacy stuff ...
    save([fitFileName '_ldat.mat'],'ldat','hd'); 
    mat2csv2Dfl(ldat,[fitFileName '_ldat.csv'],0,1,hdLong);
    ldat = ldatBak;           hd = hdBak;  % urgh ugly ...

end

disp(' '); 
disp(['Grid fit for Social Evaluation Learning Task ' fitFileName ' done.'])
t = clock; tEnd = num2str(t(2:5)); 
tStart
tEnd
save([fitFileName '_done.mat'],'tStart','tEnd','parOrdGrid','grSc','grN');

return;  % whole function.




