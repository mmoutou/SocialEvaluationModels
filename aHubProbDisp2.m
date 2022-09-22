% aHubProbDisp2 : display how the likelihoods in a are learnt over trials. 
%                 For use within spm_mdp_L_xi

% toPlot = 1;
try
 aHub = MDPs.aHub{1};
catch
    ;
end
aHDim = size(aHub);
T = aHDim(end);
za = sum(aHub,1);
za = squeeze(za);
hiRHubLik = zeros(T,2);
hiRet=2;
for t=1:T
  pa = squeeze( aHub(hiRet,:,t))' ./ za(:,t); 
  hiRHubLik(t,:) = pa(:);
end

aHubPFig1 = figure;

subplot1 = subplot(1,1,1,'Parent',aHubPFig1);
hold(subplot1,'on');

plot(hiRHubLik);
title('Learning of Lo-Hi attribution -> Positive outcome')
legend('Lo->Pos','Hi->Pos')

