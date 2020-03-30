function F = EstimateRoP(RoP,vacup,eps,ru)
% Given the vaccination uptake and efficacy and risk of unaware, provides
% the objective value for the precived R0 valuev under the true probability
% of infection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RoP - perceived R0 value
% vacup - vaccination uptake
% eps - vaccine efficacy
% ru -risk of the unaware group

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F- The objective function

% Calculate the probability of infection for the function of the true
% probability of infection
[pinf didp]=ProbInfect( RoP,vacup,eps,1,0 );
F=ru-eps*pinf; % Determines the NE for the unaware population
end
