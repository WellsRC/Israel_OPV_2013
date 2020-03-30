function [pinf dIdp] = ProbInfect( Ro,vacup,eps,gamma,pc )
%Summary: The function calculates and returns the probability of infection
%and the marginal benefit of vaccination
%Output
% pinf - The probability of infection
% dIdp - The marginal benefit of vaccination
%Input
% Ro - The basic reproductive number
% phi - The vaccine coverage
% gamma - The proportionality constant for the perceived probabilities
% pc - The choice of the probabiloty function to use
%   pc=0 : The use of the true probability of infection  pinf=1-S_inf/(1-phi)
%   pc=1 : The use of the exponential pinf=1-exp(-gamma*Ro*(1-phi))
%   pc=2 : The use of Imax pinf=gamma*Imax
phi=vacup*eps;
if(pc==0)
    % The true probability of infection
    options2=optimset('MaxIter',10^6,'MaxFunEvals',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
    sinf= lsqnonlin(@(z)-z+log(z)/Ro+1-phi-log(1-phi)./Ro,0.1*(1-phi),0,1-phi,options2);
    % Probability of infection
    pinf=1-sinf./(1-phi);
    ds=(sinf./(1-phi)).*(Ro*(1-phi)-1)./(1-Ro*sinf);
    % Marginal benefit of vaccination
    dIdp=(-ds./(1-phi)-sinf./(1-phi)^2);
    dIdp=eps*dIdp;
    if(phi>=1-1/Ro)
        pinf=0;
        dIdp=0;
    end
elseif (pc==1)
    % The perceived probability of infection using suturated  
    % Probability of infection
    pinf=Ro.*(1-phi).^(1/gamma)./(Ro+(1-phi).^(1/gamma));
    % Marginal benefit of vaccination
    dIdp=-(1/gamma).*Ro.^2.*(1-phi).^((1/gamma)-1)./(Ro+(1-phi).^(1/gamma)).^2;
    dIdp=eps*dIdp;
elseif(pc==2)
     Ro=Ro*(1+gamma);
    % The perceived probability of infection using Imax
    % Probability of infection
    pinf=(log(1/Ro)-1+Ro*(1-phi)-log(1-phi))/Ro;
    % Marginal benefit of vaccination
    dIdp=(-1+1/(Ro*(1-phi)));
    dIdp=eps*dIdp;
    if(phi>1-1/Ro)
        pinf=0;
        dIdp=0;
    end
elseif(pc==3)
    % The perceived probability of infection using saturation
    % Probability of infection
    pinf=((log(1/Ro)-1+Ro)./Ro).*(1-phi)^(1/gamma);
    % Marginal benefit of vaccination
    dIdp=-((log(1/Ro)-1+Ro)./Ro).*(1/gamma).*(1-phi).^(1/gamma-1);
    dIdp=eps*dIdp;
end
end

