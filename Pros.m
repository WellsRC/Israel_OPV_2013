rng('shuffle');
Baselineparameters
ssize=2500;
Reffv=gamrnd(11.923324191300418,0.77/11.923324191300418,ssize,1)+1;
x=betacdf(0.369,18.7916,27.8906);
dw=betainv(x+rand(ssize,1).*(1-x),18.7916,27.8906);
pl=1./1000+(1/200-1/1000).*rand(ssize,1);
pv=(1/(3*10^6))+(1/250000-1/(3*10^6)).*rand(ssize,1);
sizeg=100;
gamma=[linspace(0.05,0.2,sizeg); linspace(0,1,sizeg); linspace(0.01,0.2,sizeg)];
fr=[1];
vacupSA=zeros(sizeg,length(Reffv),3);
kappam=zeros(sizeg,length(Reffv),3);
alpha=1;
pobj=parpool(20);
rrr=1;
for ii=1:sizeg
    for kk=1:3
        parfor jj=1:length(Reffv)  
            r=(pv(jj).*dw(jj))./(pl(jj).*0.369)*rrr;
            kappa=NEkappa2(r,Reffv(jj),0.721,eps,gamma(kk,ii),kk);
            [NE1 NE2 NE3]=SolveNE(r,r,Reffv(jj),eps,kappa,omega,1,gamma(kk,ii),alpha,kk);
            kappam(ii,jj,kk)=kappa;
            vacupSA(ii,jj,kk)=omega*(alpha*NE1);
        end
    end
end

for kk=1:3
    fnum=1;
    while exist(['721VacCoverage-MinKappa-Prosocial-rho=1-pc=' num2str(kk) '-' num2str(fnum) '-Risk=' num2str(rrr) '.txt'], 'file')
           fnum=fnum+1;
    end
    f1=fopen(['721VacCoverage-MinKappa-Prosocial-rho=1-pc=' num2str(kk) '-' num2str(fnum) '-Risk=' num2str(rrr) '.txt'],'w');
    f2=fopen(['Kappa-721VacCoverage-MinKappa-Prosocial-rho=1-pc=' num2str(kk) '-' num2str(fnum) '-Risk=' num2str(rrr) '.txt'],'w');
    
    for ii=1:sizeg  
        for jj=1:length(Reffv)
            fprintf(f1,'%32.30f ',vacupSA(ii,jj,kk));
            fprintf(f2,'%32.30f ',kappam(ii,jj,kk));
        end
        fprintf(f1,';\n');
        fprintf(f2,';\n');
    end    
    fclose('all');
end

vacupSA=zeros(sizeg,length(Reffv),3);
kappam=zeros(sizeg,length(Reffv),3);
for ii=1:sizeg
    for kk=1:3
        parfor jj=1:length(Reffv)  
            r=(pv(jj).*dw(jj))./(pl(jj).*0.369)*rrr;
            kappa=NEkappa2(r,Reffv(jj),0.79,eps,gamma(kk,ii),kk);
            [NE1 NE2 NE3]=SolveNE(r,r,Reffv(jj),eps,kappa,omega,1,gamma(kk,ii),alpha,kk);
            kappam(ii,jj,kk)=kappa;
            vacupSA(ii,jj,kk)=omega*(alpha*NE1);
        end
    end
end

for kk=1:3
    fnum=1;
    while exist(['79VacCoverage-MinKappa-Prosocial-rho=1-pc=' num2str(kk) '-' num2str(fnum) '-Risk=' num2str(rrr) '.txt'], 'file')
           fnum=fnum+1;
    end
    f1=fopen(['79VacCoverage-MinKappa-Prosocial-rho=1-pc=' num2str(kk) '-' num2str(fnum) '-Risk=' num2str(rrr) '.txt'],'w');
    f2=fopen(['Kappa-79VacCoverage-MinKappa-Prosocial-rho=1-pc=' num2str(kk) '-' num2str(fnum) '-Risk=' num2str(rrr) '.txt'],'w');
    
    for ii=1:sizeg  
        for jj=1:length(Reffv)
            fprintf(f1,'%32.30f ',vacupSA(ii,jj,kk));
            fprintf(f2,'%32.30f ',kappam(ii,jj,kk));
        end
        fprintf(f1,';\n');
        fprintf(f2,';\n');
    end    
    fclose('all');
end

delete(pobj);
clear;