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
pobj=parpool(20);
rrr=1;
    for ii=1:sizeg
        for kk=1:3
            parfor jj=1:length(Reffv)
                r=(pv(jj).*dw(jj))./(pl(jj).*0.369)*rrr;
                [NE1 NE2 NE3]=SolveNE(r,r,Reffv(jj),eps,0,omega,0,gamma(kk,ii),0,kk);
                vacupSA(ii,jj,kk)=omega*NE3;
            end
        end
    end
    delete(pobj);
    for kk=1:3
        fnum=1; 
        while exist(['721VacCoverage-Individualistic-pc=' num2str(kk) '-' num2str(fnum) '-Risk=' num2str(rrr) '.txt'], 'file')
               fnum=fnum+1;
        end
        f1=fopen(['721VacCoverage-Individualistic-pc=' num2str(kk) '-' num2str(fnum) '-Risk=' num2str(rrr) '.txt'],'w');
        for ii=1:sizeg  
            for jj=1:length(Reffv)
                fprintf(f1,'%32.30f ',vacupSA(ii,jj,kk));
            end
            fprintf(f1,';\n');
        end
        fclose('all');
    end

delete(pobj);
clear;
