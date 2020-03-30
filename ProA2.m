%% Runs the anaylsis for Figure1 in the PNAS manuscript
rng('shuffle');
Baselineparameters;
[rap rup]= Estrisk;
ssize=2500;
Reffv=gamrnd(11.923324191300418,0.77/11.923324191300418,ssize,1)+1;
x=betacdf(0.369,18.7916,27.8906);
dw=betainv(x+rand(ssize,1).*(1-x),18.7916,27.8906);
pl=1./1000+(1/200-1/1000).*rand(ssize,1);
pv=(1/(3*10^6))+(1/250000-1/(3*10^6)).*rand(ssize,1);
sizeg=100;
gamma=[linspace(0.05,0.2,sizeg); linspace(0,1,sizeg); 10.^linspace(-3.5,-1,sizeg)];
vacupSA=zeros(sizeg,length(Reffv),3);
NE1v=zeros(sizeg,length(Reffv),3);
NE2v=zeros(sizeg,length(Reffv),3);
NE3v=zeros(sizeg,length(Reffv),3);
kappam=zeros(sizeg,length(Reffv),3);

alpha=0.548;
rho=0.655/0.94;
% pobj=parpool(20);
for ii=1:sizeg
    for kk=1:1
        parfor jj=1:length(Reffv)  
            r=(pv(jj).*dw(jj))./(pl(jj).*0.369);
            ru=r;
            ra=r;
            kappa=NEkappa2(ra,Reffv(jj),0.721,eps,gamma(kk,ii),kk);
            [NE1 NE2 NE3]=SolveNE(ra,ru,Reffv(jj),eps,kappa,omega,rho,gamma(kk,ii),alpha,kk);
            kappam(ii,jj,kk)=kappa;
            vacupSA(ii,jj,kk)=omega*(alpha*rho*NE1+alpha*(1-rho)*NE2+(1-alpha)*NE3);
            NE1v(ii,jj,kk)=NE1;
            NE2v(ii,jj,kk)=NE2;
            NE3v(ii,jj,kk)=NE3;
        end
    end
end
% delete(pobj);
for kk=1:1
    fnum=1;
    while exist(['721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=' num2str(rho) '-pc=' num2str(kk) '-' num2str(fnum) '.txt'], 'file')
           fnum=fnum+1;
    end
    f1=fopen(['721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=' num2str(rho) '-pc=' num2str(kk) '-' num2str(fnum) '.txt'],'w');
    f2=fopen(['Kappa-721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=' num2str(rho) '-pc=' num2str(kk) '-' num2str(fnum) '.txt'],'w');
    f3=fopen(['ProStrat-721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=' num2str(rho) '-pc=' num2str(kk) '-' num2str(fnum) '.txt'],'w');
    f4=fopen(['IndStrat-721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=' num2str(rho) '-pc=' num2str(kk) '-' num2str(fnum) '.txt'],'w');
    f5=fopen(['UnawareStrat-721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=' num2str(rho) '-pc=' num2str(kk) '-' num2str(fnum) '.txt'],'w');
    
    for ii=1:sizeg  
        for jj=1:length(Reffv)
            fprintf(f1,'%32.30f ',vacupSA(ii,jj,kk));
            fprintf(f2,'%32.30f ',kappam(ii,jj,kk));
            fprintf(f3,'%32.30f ',NE1v(ii,jj,kk));
            fprintf(f4,'%32.30f ',NE2v(ii,jj,kk));
            fprintf(f5,'%32.30f ',NE3v(ii,jj,kk));
        end
        fprintf(f1,';\n');
        fprintf(f2,';\n');
        fprintf(f3,';\n');
        fprintf(f4,';\n');
        fprintf(f5,';\n');
    end    
    fclose('all');
end
clear;
ProA;