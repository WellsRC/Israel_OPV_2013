function F=NEAware(p1,p2,p3,ra,ru,Ro,eps,kappa,omega,rho,gamma,alpha,player,pc)
%Calculates and outputs the payoff to vaccinate minus the payoff not to
%vaccinate for the specified player
%Output
% F - The payoff to vaccinate minus the payoff not to
%vaccinate for the specified player
%Input
% p1 - The vaccination strategy for player 1 (i.e. the aware prosocial group)
% p2 - The vaccination strategy for player 2 (i.e. the aware selfish group)
% p3 - The vaccination strategy for player 3 (i.e. the unaware selfish group)
% ra - The perceived relative risk for the aware group
% ru - The perceived relative risk for the unaware group
% Ro - The basic reproductive number 
% eps - The vaccine efficacy
% kappa - The relative value of someone else's healthc comapred to your own
% omega - The fraction of the population elgible for OPV vaccination 
% rho - The fraction of the aware group that is acting prosoically
% gamma - The constant required for calculating the perceived probability of infection 
% alpha - The fraction of the population aware of the prosocial nature of the campaign
% player - The player to calcualte the differences in payoffs
% pc - The choice for the probability of infection
%   pc=0 : The use of the true probability of infection  pinf=1-S_inf/(1-phi)
%   pc=1 : The use of the exponential pinf=1-exp(-gamma*Ro*(1-phi))
%   pc=2 : The use of Imax pinf=gamma*Imax

%Calculating the vaccine coverage
phi=omega*eps*(alpha.*rho.*p1+alpha.*(1-rho).*p2+(1-alpha)*p3);
vacup=omega*(alpha.*rho.*p1+alpha.*(1-rho).*p2+(1-alpha)*p3);
%Calcualte the probability of infection and marginal benefit
[pinf dIdp] = ProbInfect( Ro,vacup,eps,gamma,pc );
if(player==1)
    %Difference in the payoff to vaccinate and the payoff not to vaccinate
    %for the aware prosocial player
    F=-ra-kappa*(1-phi)*dIdp-(1-eps)*pinf+pinf*(1+kappa*Ro*(1-phi)); 
elseif(player==2)
    %Difference in the payoff to vaccinate and the payoff not to vaccinate
    %for the aware prosocial player
    F=-ra-(1-eps)*pinf+pinf;    
else
    %Difference in the payoff to vaccinate and the payof not to vaccinate
    %for the unaware selfish player
    F=-ru-(1-eps)*pinf+pinf;
end
end

