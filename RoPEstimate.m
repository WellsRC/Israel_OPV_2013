function y = RoPEstimate(ru,vacup,eps)
% runs the optimizatino algorithm to produce the estimated perceived R0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ru -risk of the unaware group
% vacup - vaccination uptake
% eps - vaccine efficacy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y- The perceived basic reproductive number R0p

% specifications for the algorithm
options=optimset('MaxFunEvals',10^6,'MaxIter',10^6,'TolFun',10^(-8),'TolX',10^(-8),'display','off');
% Runs the optimization algorithm to determine R0p
y=lsqnonlin(@(x)EstimateRoP(x,vacup,eps,ru),8,1/(1-vacup*eps),inf,options);
end

