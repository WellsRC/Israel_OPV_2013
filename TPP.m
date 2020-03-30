clear;
rng('shuffle');
Baselineparameters
ssize=2500;
x=betacdf(0.369,18.7916,27.8906);
dw=betainv(x+rand(ssize,1).*(1-x),18.7916,27.8906);
pl=1./1000+(1/200-1/1000).*rand(ssize,1);
pv=(1/(3*10^6))+(1/250000-1/(3*10^6)).*rand(ssize,1);
sizeg=50;
Reffv=linspace(1.8,2.2,sizeg);
vacupSA=zeros(sizeg,length(dw));
kappam=zeros(sizeg,length(dw));
alpha=1;
pobj=parpool(20);
for ii=1:sizeg
    parfor jj=1:length(dw)
            r=(pv(jj).*dw(jj))./(pl(jj).*0.369);
            kappa=NEkappa2(r,Reffv(ii),0.721,eps,0,0);
            [NE1 NE2 NE3]=SolveNE(r,r,Reffv(ii),eps,kappa,omega,1,0,alpha,0);
            vacupSA(ii,jj)=omega*(alpha.*NE1);
            kappam(ii,jj)=kappa;
    end
end
delete(pobj);
f2=fopen(['VacCoverage-TrueProbInfection-Prosocial.txt'],'w');
f1=fopen(['VacCoverage-Kappa-TrueProbInfection-Prosocial.txt'],'w');

    for ii=1:sizeg   
        for jj=1:length(dw)
            fprintf(f2,'%32.30f ',vacupSA(ii,jj));
            fprintf(f1,'%32.30f ',kappam(ii,jj));
        end        
        fprintf(f2,';\n');
        fprintf(f1,';\n');
    end
fclose('all');
clear;
