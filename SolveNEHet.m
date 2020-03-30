function [ TNE1 TNE2 TNE3 TNE4] = SolveNEHet(ra,ru,Ro,eps,kappaa,kappau,omega,rho,rhou,gamma,alpha,pc)
%Summary: The algorithm used to calcualte the Nash Equilibrium for the game
%between prosocial and self-interest individuals
%Outpout
% NE1 - The Nash equilibrium strategy for player 1 (i.e. the aware prosocial
% individual)
% NE2 - The Nash equilibrium strategy for player 2 (i.e. the aware self-interest
% individual)
% NE3 - The Nash equilibrium strategy for player 3 (i.e. the unaware self-interest
% individual)
% NE4 - The Nash equilibrium strategy for player 4 (i.e. the unaware
% prosocial
% individual)
% Input
% ra - The perceived relative risk for the aware group
% ru - The perceived relative risk for the unaware group
% Ro - The basic reproductive number 
% eps - The vaccine efficacy
% kappaa - The relative value of someone else's healthc comapred to your own for aware group
% kappau - The relative value of someone else's healthc comapred to your own for previously unaware group
% omega - The fraction of the population elgible for OPV vaccination 
% rho- The fraction of the aware population acting prosocially
% rhou- The fraction of the unaware population now acting prosocially
% gamma - The constant required for calculating the perceived probability of infection 
% alpha - The fraction of the population aware of the prosocial nature of the campaign
% player - The player to calcualte the differences in payoffs
% pc - The choice for the probability of infection
%   pc=0 : The use of the true probability of infection  pinf=1-S_inf/(1-phi)
%   pc=1 : The use of the exponential pinf=1-exp(-gamma*Ro*(1-phi))
%   pc=2 : The use of Imax pinf=gamma*Imax

% If the relative risk is the same among the aware and unaware then the
% self-interest group will act the same for the ware an unaware. Therefore,
% the fraction of aware prosocial becomes alpha*rho and the
% self-interest is alpha*(1-rho)+(1-alpha). Thus setting alpha to
% alpha*rho and rho to one we reduce the game to two players
if((ra==ru)&&(kappaa==kappau))
   ps=rho*alpha+(1-alpha)*rhou;
   alpha=1;
   rho=ps;
elseif((ra==ru)&&(kappaa==0))
    rho=0;
    alpha=0;
elseif((ra==ru)&&(kappau==0))
    rhou=0;
    alpha=1;
elseif(ra==ru)
    alpha=alpha*rho;
    rho=1;
elseif(kappaa==0)
    rho=0;
elseif(kappau==0)
    rhou=0;
end
options=optimset('MaxFunEvals',10^6,'MaxIter',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
xs=(1-1/Ro)/eps*0.95;
%Determine who to solve first
player=1;
%Solve the strategy for a purely aware prosocial population
NE1=lsqnonlin(@(y)NEAwareHet(y,0,0,0,ra,ru,Ro,eps,kappaa,kappau,omega,1,0,gamma,1,player,pc),xs,[0],[1],options);
player=2;
%Solve the strategy for a purely aware self-interst population
NE2=lsqnonlin(@(y)NEAwareHet(0,y,0,0,ra,ru,Ro,eps,kappaa,kappau,omega,0,0,gamma,1,player,pc),xs,[0],[1],options);
player=3;
%Solve the strategy for a purely unaware self-interst population
NE3=lsqnonlin(@(y)NEAwareHet(0,0,y,0,ra,ru,Ro,eps,kappaa,kappau,omega,0,0,gamma,0,player,pc),xs,[0],[1],options);
player=4;
%Solve the strategy for a purely unaware prosocial population
NE4=lsqnonlin(@(y)NEAwareHet(0,0,0,y,ra,ru,Ro,eps,kappaa,kappau,omega,0,1,gamma,0,player,pc),xs,[0],[1],options);
%If the strategy for the prosocial individual is greater than that of the
%self interest then the prosocial individuals desire higher vaccine
%coverage which allows the self-interest ot free ride. Otherwise, the
%self-interst population desires the higher coverage and prosocial people
%will free ride on the self interest.
%% Analytical soluation
% If the strategy for the prosocial individual is greater than that of the
% self interest then the prosocial individuals desire higher vaccine
% coverage which allows the self-interest ot free ride. Otherwise, the
% self-interst population desires the higher coverage and prosocial people
% will free ride on the self interest.
if(max([NE1 NE2 NE3 NE4])==NE1)
    TNE1=min(NE1./(alpha*rho),1);
    if(max([NE2 NE3 NE4])==NE2)
        TNE2=min(max((NE2-rho.*alpha*TNE1)./(alpha*(1-rho)),0),1);  
        if(max([NE3 NE4])==NE3)
            TNE3=min(max((NE3-rho.*alpha*TNE1-alpha*(1-rho)*TNE2)./((1-rhou)*(1-alpha)),0),1);   
            TNE4=min(max((NE4-rho.*alpha*TNE1-alpha*(1-rho)*TNE2-TNE3*(1-alpha)*(1-rhou))./(rhou*(1-alpha)),0),1);
        else
            TNE4=min(max((NE4-rho.*alpha*TNE1-alpha*(1-rho)*TNE2)./(rhou*(1-alpha)),0),1);
            TNE3=min(max((NE3-rho.*alpha*TNE1-alpha*(1-rho)*TNE2-TNE4*rhou*(1-alpha))./((1-rhou)*(1-alpha)),0),1); 
        end
    elseif(max([NE2 NE3 NE4])==NE3)            
        TNE3=min(max((NE3-rho.*alpha*TNE1)./((1-rhou)*(1-alpha)),0),1);
        if(max([NE2 NE4])==NE2)
            TNE2=min(max((NE2-rho.*alpha*TNE1-(1-alpha)*(1-rhou)*TNE3)./(alpha*(1-rho)),0),1);  
            TNE4=min(max((NE4-rho.*alpha*TNE1-alpha*(1-rho)*TNE2-(1-alpha)*(1-rhou)*TNE3)./(rhou*(1-alpha)),0),1);
        else
            TNE4=min(max((NE4-rho.*alpha*TNE1-(1-alpha)*(1-rhou)*TNE3)./(rhou*(1-alpha)),0),1);   
            TNE2=min(max((NE2-rho.*alpha*TNE1-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*(1-rho)),0),1); 
        end
    else
        TNE4=min(max((NE4-rho.*alpha*TNE1)./(rhou*(1-alpha)),0),1); 
        if(max([NE2 NE3])==NE2)
            TNE2=min(max((NE2-rho.*alpha*TNE1-rhou*(1-alpha)*TNE4)./(alpha*(1-rho)),0),1);  
            TNE3=min(max((NE3-rho.*alpha*TNE1-alpha*(1-rho)*TNE2-TNE4*rhou*(1-alpha))./((1-rhou)*(1-alpha)),0),1); 
        else
            TNE3=min(max((NE3-rho.*alpha*TNE1-(1-alpha)*rhou*TNE4)./((1-rhou)*(1-alpha)),0),1);   
            TNE2=min(max((NE2-rho.*alpha*TNE1-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*(1-rho)),0),1); 
        end
    end
elseif(max([NE1 NE2 NE3 NE4])==NE2)
    TNE2=min(NE2./(alpha*(1-rho)),1);
    if(max([NE1 NE3 NE4])==NE1)
        TNE1=min(max((NE1-(1-rho).*alpha*TNE2)./(alpha*rho),0),1);  
        if(max([NE3 NE4])==NE3)
            TNE3=min(max((NE3-rho.*alpha*TNE1-alpha*(1-rho)*TNE2)./((1-rhou)*(1-alpha)),0),1);   
            TNE4=min(max((NE4-rho.*alpha*TNE1-alpha*(1-rho)*TNE2-TNE3*(1-alpha)*(1-rhou))./(rhou*(1-alpha)),0),1);
        else
            TNE4=min(max((NE4-rho.*alpha*TNE1-alpha*(1-rho)*TNE2)./(rhou*(1-alpha)),0),1);
            TNE3=min(max((NE3-rho.*alpha*TNE1-alpha*(1-rho)*TNE2-TNE4*rhou*(1-alpha))./((1-rhou)*(1-alpha)),0),1); 
        end
    elseif(max([NE1 NE3 NE4])==NE3)            
        TNE3=min(max((NE3-(1-rho).*alpha*TNE2)./((1-rhou)*(1-alpha)),0),1);
        if(max([NE1 NE4])==NE1)
            TNE1=min(max((NE1-(1-rho).*alpha*TNE2-(1-alpha)*(1-rhou)*TNE3)./(alpha*rho),0),1);  
            TNE4=min(max((NE4-(1-rho).*alpha*TNE2-alpha*(1-rho)*TNE2-(1-alpha)*(1-rhou)*TNE3)./(rhou*(1-alpha)),0),1);
        else
            TNE4=min(max((NE4-(1-rho).*alpha*TNE2-(1-alpha)*(1-rhou)*TNE3)./(rhou*(1-alpha)),0),1);   
            TNE1=min(max((NE1-(1-rho).*alpha*TNE2-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*rho),0),1); 
        end
    else
        TNE4=min(max((NE4-(1-rho).*alpha*TNE2)./(rhou*(1-alpha)),0),1); 
        if(max([NE1 NE3])==NE1)
            TNE1=min(max((NE1-(1-rho).*alpha*TNE2-rhou*(1-alpha)*TNE4)./(alpha*rho),0),1);  
            TNE3=min(max((NE3-rho.*alpha*TNE1-alpha*(1-rho)*TNE2-TNE4*rhou*(1-alpha))./((1-rhou)*(1-alpha)),0),1); 
        else
            TNE3=min(max((NE3-(1-rho).*alpha*TNE2-(1-alpha)*rhou*TNE4)./((1-rhou)*(1-alpha)),0),1);   
            TNE1=min(max((NE1-(1-rho).*alpha*TNE2-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*rho),0),1); 
        end
    end
elseif(max([NE1 NE2 NE3 NE4])==NE3)
    TNE3=min(NE3./((1-rhou)*(1-alpha)),1);
    if(max([NE1 NE2 NE4])==NE1)
        TNE1=min(max((NE1-(1-rhou)*(1-alpha)*TNE3)./(alpha*rho),0),1);  
         if(max([NE2 NE4])==NE2)
            TNE2=min(max((NE2-rho.*alpha*TNE1-(1-alpha)*(1-rhou)*TNE3)./(alpha*(1-rho)),0),1);  
            TNE4=min(max((NE4-rho.*alpha*TNE1-alpha*(1-rho)*TNE2-(1-alpha)*(1-rhou)*TNE3)./(rhou*(1-alpha)),0),1);
        else
            TNE4=min(max((NE4-rho.*alpha*TNE1-(1-alpha)*(1-rhou)*TNE3)./(rhou*(1-alpha)),0),1);   
            TNE2=min(max((NE2-rho.*alpha*TNE1-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*(1-rho)),0),1); 
        end
    elseif(max([NE1 NE2 NE4])==NE2)            
        TNE2=min(max((NE2-(1-rhou)*(1-alpha)*TNE3)./(alpha*(1-rho)),0),1); 
        if(max([NE1 NE4])==NE1)
            TNE1=min(max((NE1-(1-rho).*alpha*TNE2-(1-alpha)*(1-rhou)*TNE3)./(alpha*rho),0),1);  
            TNE4=min(max((NE4-(1-rho).*alpha*TNE2-alpha*(1-rho)*TNE2-(1-alpha)*(1-rhou)*TNE3)./(rhou*(1-alpha)),0),1);
        else
            TNE4=min(max((NE4-(1-rho).*alpha*TNE2-(1-alpha)*(1-rhou)*TNE3)./(rhou*(1-alpha)),0),1);   
            TNE1=min(max((NE1-(1-rho).*alpha*TNE2-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*rho),0),1); 
        end
    else
        TNE4=min(max((NE4-(1-rhou)*(1-alpha)*TNE3)./((1-alpha)*rhou),0),1);
        if(max([NE1 NE2])==NE1)
            TNE1=min(max((NE1-(1-alpha)*TNE3-rhou*(1-alpha)*TNE4)./(alpha*rho),0),1);  
            TNE2=min(max((NE2-rho.*alpha*TNE1-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*(1-rho)),0),1); 
        else
            TNE2=min(max((NE2-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*(1-rho)),0),1); 
            TNE1=min(max((NE1-(1-rho).*alpha*TNE2-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*rho),0),1); 
        end
    end
else
    TNE4=min(NE4./(rhou*(1-alpha)),1);   
    if(max([NE1 NE2 NE3])==NE1)
        TNE1=min(max((NE1-(rhou*(1-alpha))*TNE4)./(alpha*rho),0),1);  
         if(max([NE2 NE3])==NE2)
            TNE2=min(max((NE2-rho.*alpha*TNE1-rhou*(1-alpha)*TNE4)./(alpha*(1-rho)),0),1);  
            TNE3=min(max((NE3-rho.*alpha*TNE1-alpha*(1-rho)*TNE2-TNE4*rhou*(1-alpha))./((1-rhou)*(1-alpha)),0),1); 
        else
            TNE3=min(max((NE3-rho.*alpha*TNE1-(1-alpha)*rhou*TNE4)./((1-rhou)*(1-alpha)),0),1);   
            TNE2=min(max((NE2-rho.*alpha*TNE1-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*(1-rho)),0),1); 
        end
    elseif(max([NE1 NE2 NE3])==NE2)            
        TNE2=min(max((NE2-(rhou*(1-alpha))*TNE4)./(alpha*(1-rho)),0),1); 
        if(max([NE1 NE3])==NE1)
            TNE1=min(max((NE1-(1-rho).*alpha*TNE2-rhou*(1-alpha)*TNE4)./(alpha*rho),0),1);  
            TNE3=min(max((NE3-rho.*alpha*TNE1-alpha*(1-rho)*TNE2-TNE4*rhou*(1-alpha))./((1-rhou)*(1-alpha)),0),1); 
        else
            TNE3=min(max((NE3-(1-rho).*alpha*TNE2-(1-alpha)*rhou*TNE4)./((1-rhou)*(1-alpha)),0),1);   
            TNE1=min(max((NE1-(1-rho).*alpha*TNE2-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*rho),0),1); 
        end
    else
        TNE3=min(max((NE3-(rhou*(1-alpha))*TNE4)./((1-alpha)*(1-rhou)),0),1);
        if(max([NE1 NE2])==NE1)
            TNE1=min(max((NE1-(1-alpha)*TNE3-rhou*(1-alpha)*TNE4)./(alpha*rho),0),1);  
            TNE2=min(max((NE2-rho.*alpha*TNE1-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*(1-rho)),0),1); 
        else
            TNE2=min(max((NE2-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*(1-rho)),0),1); 
            TNE1=min(max((NE1-(1-rho).*alpha*TNE2-(1-alpha)*TNE3-(rhou*(1-alpha))*TNE4)./(alpha*rho),0),1); 
        end
    end
end

%If the risk for the aware and unaware are the smae then the selfinterest
%groups will be have the same. Sicne we solved for group 3 (i.e. alpha adjusted to alpha*rho then (1-rho) adjusted to 0), we set NE2=NE3.
if((ra==ru)&&(kappaa==kappau))
   TNE3=TNE2;
   TNE4=TNE1;
elseif((ra==ru)&&(kappaa==0))
    TNE2=TNE3;
    TNE1=TNE3;
elseif((ra==ru)&&(kappau==0))
    TNE3=TNE2;
    TNE4=TNE2;
elseif(ra==ru)
    TNE2=TNE3;
elseif(kappaa==0)
    TNE1=TNE2;
elseif(kappau==0)
    TNE4=TNE3;
end

end

