function F= Estimategamma( gamma,RoP,vacup,eps,pc,ru )
% Given the vaccination uptake and efficacy and Ro, provides
% the objective value for gamma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gamma - the values to evalaute for the objective function
% RoP - perceived R0 value
% vacup - vaccination uptake
% eps - vaccine efficacy
% pc - the proba infectino functino to use
% ru -risk of the unaware group

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F- The objective function

% Calculate the probability of infection for the function of the true
% probability of infection
[pinf didp]=ProbInfect( RoP,vacup,eps,gamma,pc );
% Minimize the NE
F=(ru)-(eps*pinf);
end

