function [rho,minavgkappa,kappam] = TableValues(ra,ru,RoP,vacup,eps,omega,alpha,gamma,pc )
%TABLEVALUES returns the estimated values required for the table in the
%results
%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%
% ra - The perceived relative risk for the aware group
% ru - The perceived relative risk for the unaware group
% RoP - The basic reproductive number 
% Vacup - vaccination uptake
% eps - The vaccine efficacy
% omega - The fraction of the population elgible for OPV vaccination 
% alpha - The fraction of the population aware of the prosocial nature of the campaign
% gamma - The constant required for calculating the perceived probability of infection 
% pc - The choice for the probability of infection
%   pc=0 : The use of the true probability of infection  pinf=1-S_inf/(1-phi)
%   pc=1 : The use of the exponential pinf=1-exp(-gamma*Ro*(1-phi))
%   pc=2 : The use of Imax pinf=gamma*Imax

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rho - proportion proscoail
% minavgkappa - minimum kappa value
% kappam - alternative estimate of kappa 

%options=optimset('MaxIter',10^6,'MaxFunEvals',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off','Algorithm','interior-point');

% determine the minimum value of kappa to be returned
kappa=NEkappa2(ra,RoP,vacup,eps,gamma,pc); % calculates the minimum kappa valuer
rho=vacup/(omega*(1.224-0.224*alpha)); % determine the proprtion who were prosocial
minavgkappa=kappa;  % sets minimum kappa value

% Determine the alternative value of kappa 
if(pc>0)
    [pinf dIdp] = ProbInfect( RoP,vacup,eps,gamma,pc );
else
    [pinf dIdp] = ProbInfect( RoP,vacup,eps,0,pc );
end
a=-dIdp*(1-eps*vacup);
b=(-ra-(1-eps)*pinf);
c=dIdp*(1-omega);
kappam=(-b+sqrt(b^2-4*a*c))./(2*a); % Estimated kappa value

end

