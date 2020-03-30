function kappa= NEkappa2(ra,Ro,vacup,eps,gamma,pc)

%Determines the value of kappa for a given vaccination uptake

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ra- risk for the aware population
% Ro - the basic reproductive number
% vacup - the vaccination uptake
% eps - the vaccine efficacy
% gamma - the risk of infection for the probability of infection
% omega - the percentage of the populatino elgible for vaccination
% pc - the probability of infection that is to be used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kappa - the levels of prosocial behaviour

% Adjustments for herd immunity functions
phi=vacup*eps; % Determine the coverage i.e. uptake x efficacy
if((pc==0)||(pc==2))
    if(phi>1-1/(Ro*(1+gamma)))
        phi=(1-1/(Ro*(1+gamma)))*0.9999;
        vacup=phi/eps;
    end
end    

% Calcualte the probability of infection and marginal benefit
[pinf,dIdp] = ProbInfect( Ro,vacup,eps,gamma,pc );
% the estimateed value of kappa
kappa=(ra-pinf*eps)/((1-phi)*(-dIdp+Ro*pinf));

% ensure that kappa is in zeor and one
if((kappa<0)||(kappa==Inf))
    kappa=0;
end
if(kappa>1)
    kappa=1;
end
end

