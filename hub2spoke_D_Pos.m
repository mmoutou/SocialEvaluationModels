function posRep_D = hub2spoke_D_Pos(hub_a1, hub_d, allP )
% function posRep_D = hub2spoke_D_Pos(hub_a1,hub_d, allP )
% Provides normalised prob. D : posRep_D{1} = upmf' / sum(upmf);  
%          Estimate the D map for reporting expected positive outcome, based on the
% effective counts in the a map of the hub. That is, we assume that the cells of
% the a map summarise all effective observations that the pt considers for
% this particular interaction.
%      hub_a1 is the appropriate trial mdpHub(trial).a{1}

  resNRepo = allP.resNRepo;
  noiseFloor = allP.noiseFloor;

  % the uncertainty in the state does not contribute to the
  % final uncertainty !!! :
  ph = hub_d{1}/ sum(hub_d{1}) ;
  pAttrLH = ph ;
  
  fa = squeeze(hub_a1(2,:));         % favourable, or positive, outcomes
  ua = squeeze(hub_a1(1,:));         % unfavourable, or negative, outcomes
  ma = fa ./ (fa + ua);            % means of a
  na = (ua+fa);                    % notional evidence for each state
  va = (fa .* ua) ./ ((na.^2) .* (na+1));   % variances within a 
  
  M  = sum(sum( pAttrLH .* ma' ));    % Grand mean
  S2 = (M - ma).^2; 
  Va = S2 + va;   % variances: sum of indiv. variances plus means deviates squared
  V  = Va' .* pAttrLH ;       
  SD  = sqrt(sum(V(:)));    % overall standard dev, weighed by sample sizes.
  % The standard error of the mean is still the SD / sqrt(effective N)
  % Note that in this the effective N increase with less than 1 for each 
  % new observation !! :
  N = sum(sum(na' .* pAttrLH )); % N = length(V(:)) * sum(sum(na .* pAttrLH )); 
  SE = SD / sqrt(N);
  
  
  % Moment-matched beta distro with the mean and standard error of the 
  % mean estimate: 
  [A, B] = betaMS2ab (M, SE);
  
  cbd = [0 betacdf( (1:(resNRepo-1))/resNRepo, A, B) 1];
  upmf = (cbd(2:end) - cbd(1:resNRepo)) + noiseFloor ;     % add a bit of noise
  posRep_D{1} = upmf' / sum(upmf);                    % normalise

  % F2 are 'report of level of Intent', certain to be at their initial state:
  posRep_D{2} = [ones(resNRepo,1); 512]/(resNRepo+512);

return

