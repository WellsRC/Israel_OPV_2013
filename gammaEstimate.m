function [y] = gammaEstimate( Ro,vacup,eps,pc,ru )
% Runs the fitting algorithm to estimate gamma in the model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ro - the basic reproductive number
% vacup - the vaccination uptake
% eps - the vaccine efficacy
% gamma - the risk of infection for the probability of infection
% pc - the probability of infection that is to be used
% ra- risk for the unaware population

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y - the estimated value of gamm

options=optimset('MaxFunEvals',10^6,'MaxIter',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
if(pc~=2)
[y,f]=lsqnonlin(@(x)Estimategamma( x,Ro,vacup,eps,pc,ru ),0.1,0,inf,options);
[y2,f2]=lsqnonlin(@(x)Estimategamma( x,Ro,vacup,eps,pc,ru ),3,0,inf,options);
else
    
[y,f]=lsqnonlin(@(x)Estimategamma( x,Ro,vacup,eps,pc,ru ),0.1,-inf,inf,options);
[y2,f2]=lsqnonlin(@(x)Estimategamma( x,Ro,vacup,eps,pc,ru ),-0.1,-inf,inf,options);
end
% see which one fit best
if(f2<f)
    y=y2;
end
end

