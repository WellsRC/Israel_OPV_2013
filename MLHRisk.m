function F=MLHRisk(x)
% Used to estimate the risk of infection for the aware and unaware
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% x - the risk of infection for aware (x(1)) and unaware (x(2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% F- the negative of the log-lieklihood
F=-log((x(2)-x(1)).*(LHrisk(x(1))).*(LHrisk(x(2))));%-log((LHrisk(x(1)))); % the negative is because we want to maximize the distance between the two risk but also maximize the likelihood of each 
end

