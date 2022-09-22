function [outp,lab] = mdp2arr(mdpv, v)
% [outp,lab] = mdp2arr(mdpv, v)
% returns in array form observations if v='o', states of v='s', priorOverState for 'd',
% 
% E.g. to get a col vec of the obs. along factor 3 at time 2 for all trials, use:
% x = mdp2arr(mmdp,'o'); squeeze(x(3,2,:))'
% However if the variable is itself a cell array such as d, outp will be 
%    a cell array of same length ... 

try v; catch v = '-'; end
l = length(mdpv);
outc = {};  % temp cell array
% If observations o are required
if v == 'o'
   lab = 'observations';
   di = size(mdpv(1).o);
   outp = nan([di,l]);
   for tr = 1:l
     for ro = 1:di(1)
       outp(ro,:,tr) = mdpv(tr).o(ro,:); 
     end
   end
% if actions u are required
elseif v == 'u'
   lab = 'actions';
   di = size(mdpv(1).u);
   outp = nan([di,l]);
   for tr = 1:l
     for ro = 1:di(1)
       outp(ro,:,tr) = mdpv(tr).u(ro,:); 
     end
   end
% if states s are required
elseif v == 's'
   lab = 'trueStates';
   di = size(mdpv(1).s);
   outp = nan([di,l]);
   for tr = 1:l
     for ro = 1:di(1)
       outp(ro,:,tr) = mdpv(tr).s(ro,:); 
     end
   end
elseif v == 'd'
   lab = 'priorConcOverStates';
   dlen = length(mdpv(1).d);     % number of state factors in d
   outp = {};
   for di = 1:dlen
       ddi = length(mdpv(1).d{di});
       outp{di} = nan([ddi,l]);
       for tr = 1:l
          outp{di}(:,tr) = mdpv(tr).d{di}(:); 
       end
   end
elseif v == 'D'
   lab = 'priorPyOverStates';
   dlen = length(mdpv(1).D);     % number of state factors in D
   outp = {};
   for di = 1:dlen
       ddi = length(mdpv(1).D{di});
       outp{di} = nan([ddi,l]);
       for tr = 1:l
          outp{di}(:,tr) = mdpv(tr).D{di}(:); 
       end
   end
elseif v == 'a'
   lab = 'likelihoodConc';
   aLen = length(mdpv(1).a);     % number of state factors in a
   outp = {};
   for ai = 1:aLen
       aai = size(mdpv(1).a{ai});
       outp{ai} = [];
       for tr = 1:l  % Here we will first concatenate linearized, then reshape.
          outp{ai}= [outp{ai} mdpv(tr).a{ai}(:)]; 
       end
       outp{ai} = reshape(outp{ai},[aai,l]);
   end
else
    error(['mdp2arr not ready for v being ', v]);
end

return;