function [ra,ru]= Estrisk
% Runs the algorithm to estimate the rish for aware an unaware

%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ra- risk for aware
% ru - risk for unaware

options=optimset('MaxIter',10^6,'MaxFunEvals',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off','Algorithm','interior-point');
% Runs  the alorithm
[y]=fmincon(@(x)MLHRisk(x),[0.0003 0.0006],[],[],[],[],[0 0],[inf inf],[],options);
% Sets th estimated values
ru=y(1);
ra=y(2);
end

