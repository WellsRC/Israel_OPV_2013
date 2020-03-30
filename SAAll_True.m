ssize=5000;
colors = [255 255 191; 69 117 180; 145 191 219; 224 243 248; 254 224 144; 252 141 89; 215 48 39]./256;
Baselineparameters;
eps=0.56+0.44.*rand(ssize,1);
omega=0.9+0.08.*rand(ssize,1);
alpha=0.448+0.2.*rand(ssize,1);
x=betacdf(0.369,18.7916,27.8906);
dw=betainv(x+rand(ssize,1).*(1-x),18.7916,27.8906);
pv=1/(3*10^6)+(1/250000-1/(3*10^6)).*rand(ssize,1);
pl=1/1000+(1/200-1/1000).*rand(ssize,1);
fra=1+4.*rand(ssize,1);
V=zeros(6,3,3);
M=zeros(6,3,3);
xs=[0.075 0.377 0.692];
ys=[0.68 0.375 0.075];
wt=0.28;
ht=0.25;
vacup=0.721;
pobj=parpool(4);
RoP=zeros(ssize,6);
rho=RoP;
kappam=RoP;
minavgkappa=kappam;
ru=(pv.*dw)./(pl.*0.369);
ra=fra.*ru;
options2=optimset('MaxIter',10^8,'MaxFunEvals',10^8,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
pc=0;

    jj=1;
    for ii=1:ssize
        RoP(ii,jj) = RoPEstimate(ru(ii),0.721,eps(ii));
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(ra(ii),ru(ii),RoP(ii,jj),vacup,eps(ii),omega(ii),alpha(ii),0,pc );
    end

    %Fix eps at 0.63
    jj=2;
    parfor ii=1:ssize
        RoP(ii,jj) = RoPEstimate(ru(ii),0.721,0.63);
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(ra(ii),ru(ii),RoP(ii,jj),vacup,0.63,omega(ii),alpha(ii),0,pc );
    end
    %Fix omega at 0.93
    jj=3;
    parfor ii=1:ssize
        RoP(ii,jj) = RoPEstimate(ru(ii),0.721,eps(ii));
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(ra(ii),ru(ii),RoP(ii,jj),vacup,eps(ii),0.94,alpha(ii),0,pc );
    end
    %Fix ru at 2*10^(-4)
    jj=4;
    parfor ii=1:ssize
        RoP(ii,jj) = RoPEstimate(2*10^(-4),0.721,eps(ii));
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(2*10^(-4).*fra(ii),2*10^(-4),RoP(ii,jj),vacup,eps(ii),omega(ii),alpha(ii),0,pc );
    end
    %Fix fra at 5
    jj=5;
    parfor ii=1:ssize
        RoP(ii,jj) = RoPEstimate(ru(ii),0.721,eps(ii));
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(ru(ii).*5,ru(ii),RoP(ii,jj),vacup,eps(ii),omega(ii),alpha(ii),0,pc );
    end
    %Fix alpha
    jj=6;
    parfor ii=1:ssize
        RoP(ii,jj) = RoPEstimate(ru(ii),0.721,eps(ii));
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(ra(ii),ru(ii),RoP(ii,jj),vacup,eps(ii),omega(ii),0.548,0,pc );
    end

    for jj=1:6
        f1=fopen(['Sensitivity to infection-pc=' num2str(pc) '-Case=' num2str(jj) '.txt'],'w');
        f2=fopen(['Density of prosocials-pc=' num2str(pc) '-Case=' num2str(jj) '.txt'],'w');
        f3=fopen(['MinKappa-pc=' num2str(pc) '-Case=' num2str(jj) '.txt'],'w');
        f4=fopen(['Kappa-pc=' num2str(pc) '-Case=' num2str(jj) '.txt'],'w');
        for ii=1:ssize
            fprintf(f1,'%32.31f \n',RoP(ii,jj));
            fprintf(f2,'%32.31f \n',rho(ii,jj));
            fprintf(f3,'%32.31f \n',max([min([minavgkappa(ii,jj),1]),0]));
            fprintf(f4,'%32.31f \n',max([min([kappam(ii,jj),1]),0]));
        end
        fclose('all');
    end
delete(pobj);
