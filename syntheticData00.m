% script syntheticData00.m
% To make synthetic data for testing model-fitting and parameter recovery.
% Change base directory (folder) baseDir to use this.
% Will create if necessary, and make a subdir to write within, dataFits/fittingTests
%
% Works well enough as per tests to 1 June 23, but labelled 00 as it uses the same
% stimuli for all generated data - would be better to use stimuli generated 
% either by constrained randomisation or as seen by a set of real participants.

%% Menu-like items
codeTest = 0;  % set to 1 to activate debugging printouts etc.
synDatName = 'synRepDictCIs';   % name for synthetic dataset e.g. one with
   % ranges (confidence-interval like) deduced from a repeated dictator study.
synTotN = 1000;  
if codeTest; synTotN =3; end

numIDfiller = 660000; 

%  Set High and Lo end of range Parameter vector in native space. Set them really close
%  together, differing by 'little', if you want that param effectively constant
little = 1e-5;  % number which should be psychologically negligible 
%         1          2        3      4         5        6       7        8        9    10    11     12      
%      posiSelf  posiOther dEvSelf dEvOther aEvSelf aEvOther alphaPrec genLR    repLR wp0   wAttr  mem      
parHi =  [0.9,     0.8,     40,    40,     2+little,2+little,  10, 0.1+little,  0.5,  1.26,  10,   0.999+little ]; 
parLo =  [0.4,     0.3,     2,     2,      2-little,2-little,   5, 0.1-little,  0.1, -1.1,   0.05, 0.999-little ]; 
% parHi =  [0.9,     0.8,     507.8,  507.8,   70.53,   70.53,  10.78,    0.2,   0.3, 1.26,  10,   0.999+little ]; 
% parLo =  [0.4,     0.3,     2.12,   2.12,    0.042,   0.042,  0.690,    0.01, 0.05, -1.1,  0.05, 0.999-little ]; 
% Same in transformed space:
pTrHi = nat2tr_SEL_i(parHi,1);  % map to transformed space. Transforms all the params, not just for each block.
pTrLo = nat2tr_SEL_i(parLo,1);  % rem last entry means we are working with vectors not structures for now.
pTrM  = (pTrHi + pTrLo)/2;
pTrSD = abs(pTrHi - pTrLo)/4;   % NB ranges above should corresp. to +/- 2SD

%% The path where dataFits/fittingTests resides here.
%  Make note of where we are and change to it.
baseDir = {};
baseDir{1} = 'NewUserBasePath';   % Users can add their path here. 
baseDir{2} = '/home/mmoutou/googledirs/MM/SocEvalModels_(not_git)_rough_work/';
baseDir{3} = 'C:\Users\mmpsy\Nextcloud\MM\googledirs\SocEvalModels_(not_git)_rough_work\';
% not sure if this exists on Ubuntu PC :
baseDir{4} = '/home/michael/googledirs/MM/SocEvalModels_(not_git)_rough_work/';
cwd = pwd;
baseD = {};
for dirN = 1:length(baseDir)
    if isempty(baseD)
      try
        cd(baseDir{dirN});
        baseD=  baseDir{dirN} ; 
      catch
        baseD = {};
        warning([baseDir{dirN} ' N/A ...']);
      end
    end
end

%% A bit of housework
hd=  {'ptID','posiSelf','posiOther','dEvSelf','dEvOther','aEvSelf','aEvOther','alphaPrec','genLR','repLR','wp0','wAttr','mem','LL'};
hdPar=  {'posiSelf','posiOther','dEvSelf','dEvOther','aEvSelf','aEvOther','alphaPrec','genLR','repLR','wp0','wAttr','mem'};
hd2=  {'ptID','posiSelf','posiOther','dEvSelf','dEvOther','aEvSelf','aEvOther','alphaPrec','genLR','repLR','wp0','wAttr','mem'};
repoFun= 1;   % 1 for simple softmax-like;
load('noLearnSEL8bl.mat');  % load modStruc8, mdp8, inp8, resp8,pPosGen,par8hd,selfpHd ,
    % several of which will be used either in generative model or as templates for outputs.
details = 1;    % produce detailed output, but don't plot
% fields:         posiSelf posiOther dEvSelf dEvOther aEvSelf aEvOther alphaPrec genLR repLR wp0   wAttr  mem    
modStruc8.priPar=[[1.01,  1.01,     1.01,    1.01,  1.01,    1.01,     1.01,   1.01, 1.01, 10,    1.2,  1.01]; ...  % A
                    [1.01,  1.01,      2,       2,      2,       2,        2,    1.01, 1.01, 10,  5.8,  1.01]; ...  % B
                    [ 0,     0,        0,       0,      0,       0,        0,     0,   0,   -46,   0,    0   ]; ...  % lo
                    [ 1,     1,        100,    100,    100,     100,     100,     1,   1,    46   100,   1   ]];     % hi
%                                    <-  max             at 1                >            <SD ~10> <max at 4>
modStruc8.repoFun = repoFun;

%% Sample the parameters in transf. space
%  In this rather frugal example, set dEvOther = dEvSelf ; aEvOther = aEvSelf; 
%   
% params that will have same values:
yokePar  = [[3 4];  [5 6]];          % 3. determines 4.,  5. determines 6.
totParN  = length(parHi);
genTrPar = nan(synTotN,totParN);     % to store transf. space generative params.
genNatPar= genTrPar;                 % to store native space gen params.
genNatParLL = nan(synTotN,totParN+2); % to store transf. space generative params w. an ID and the LL.

for ptN = 1:synTotN
    for parN = 1:totParN
        found = find(yokePar(:,2)==parN);
        if found
           genTrPar(ptN,parN) = genTrPar(ptN,yokePar(found,1));
        else
           genTrPar(ptN,parN) = normrnd(pTrM(parN),pTrSD(parN));
        end
    end
    genNatPar(ptN,:) = tr2nat_SEL_i(genTrPar(ptN,:));  % does not need other arguments :o
    genNatParLL(ptN,1) = numIDfiller + ptN;
    genNatParLL(ptN,2:(end-1)) = genNatPar(ptN,:);
end
         

%% main loop to make synthetic data
selD = {};
for ptN = 1:synTotN
    disp(['Now on ' num2str(ptN) ' out of ' num2str(synTotN)]);
    pTr8 = genTrPar(ptN,:); 
    Inp = inp8;   % Ideally should come from actual experiment, or be pseudorandom.
    Resp = [];    % v. important to elicit generative mode in next line.
    [syFitM, syPar, syMDP8] = sel8bl03( pTr8, Inp, Resp, modStruc8, repoFun, details);
    Resp = resp8;    % to use as template.
    for blN=1:length(Inp)
       for trN=1:length(Inp{blN}.posRep)
           Resp{blN}.posRep{trN} = syMDP8{blN}.Resp.posRep{trN}; 
       end
    end 
    selD{ptN}.trPar = pTr8;
    selD{ptN}.Inp = Inp;
    selD{ptN}.Resp = Resp;
    selD{ptN}.fitM = syFitM;
    selD{ptN}.MDP8 = syMDP8; 
    selD{ptN}.numID = ptN + numIDfiller;

    genNatParLL(ptN,end) = selD{ptN}.fitM.sLL; 
  
end

%% Storage directory -------------------------------------
if exist( 'dataFits/fittingTests','dir')
    cd('dataFits/fittingTests');
else
    warning('dataFits/fittingTests does not exist -- working here:');
    pwd
end
mkdir(pwd,synDatName);
cd(synDatName);

%% Save the model and sampled parameters in mat and csv form 
save('modStruc8a.mat','modStruc8');  % will be looked for by e.g. SEL2_Grid04a_xiii
%  Params in native space, flanked by IDs and LLs
nameFrag = [num2str(synTotN) '_' synDatName];
mat2csv2Dfl(genNatParLL,['genNatParLL_' nameFrag '.csv'],0,1,hd); 
%  Same in transf. space, to use in re-fitting / recovery. HAS TO HAVE 'hd', not 'hd2' etc.
ldat = [genNatParLL(:,1) genTrPar genNatParLL(:,end)]; 
save(['ldat_gen_' nameFrag '.mat'], 'ldat', 'hd');
mat2csv2Dfl(ldat,['ldat_gen_' nameFrag '.csv'],0,1,hd); 

%% Save the synthetic data and wrap up ------------------
save(['dat' nameFrag '.mat'],'selD','modStruc8','ldat','hd','genNatParLL','hdPar');

disp(['Synthetic data saved in ' pwd]);
disp(['for run' nameFrag]);
disp('Generative params genNatPar and genTrPar in respective mat and csv files.');
disp('In this version, all synthetic pts have same input');
disp('explore with e.g. [fitMsy2, Par8sy2, MDP8sy2] = sel8bl03( selD{2}.trPar, selD{2}.Inp, selD{2}.Resp,');
disp('                                                          modStruc8, repoFun, details);');


cd( cwd ); 



% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  end of file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


