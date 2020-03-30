ssize=5000;
colors = [255 255 191; 69 117 180; 145 191 219; 224 243 248; 254 224 144; 252 141 89; 215 48 39]./256;
Baselineparameters;
eps=0.56+0.44.*rand(ssize,1);
omega=0.9+0.08.*rand(ssize,1);
RoP=gamrnd(11.923324191300418,0.77/11.923324191300418,ssize,1)+1;
alpha=0.448+0.2.*rand(ssize,1);
x=betacdf(0.369,18.7916,27.8906);
dw=betainv(x+rand(ssize,1).*(1-x),18.7916,27.8906);
pv=1/(3*10^6)+(1/250000-1/(3*10^6)).*rand(ssize,1);
pl=1/1000+(1/200-1/1000).*rand(ssize,1);
fra=1+2.*rand(ssize,1);
V=zeros(6,3,3);
M=zeros(6,3,3);
xs=[0.075 0.377 0.692];
ys=[0.68 0.375 0.075];
wt=0.28;
ht=0.25;
vacup=0.721;

gamma=zeros(ssize,7);
rho=gamma;
kappam=gamma;
minavgkappa=kappam;
ru=(pv.*dw)./(pl.*0.369);
ra=fra.*ru;
options2=optimset('MaxIter',10^8,'MaxFunEvals',10^8,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
 for pc=1:3

    jj=1;
   parfor ii=1:ssize
        gamma(ii,jj) = gammaEstimate( RoP(ii),0.721,eps(ii),pc,ru(ii) );
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(ra(ii),ru(ii),RoP(ii),vacup,eps(ii),omega(ii),alpha(ii),gamma(ii,jj),pc );
    end
    %Fix RoP
    jj=2;
    parfor ii=1:ssize
        gamma(ii,jj) = gammaEstimate(2.25,0.721,eps(ii),pc,ru(ii) );
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(ra(ii),ru(ii),2.25,vacup,eps(ii),omega(ii),alpha(ii),gamma(ii,jj),pc );
    end

    %Fix eps at 0.63
    jj=3;
    parfor ii=1:ssize
        gamma(ii,jj) = gammaEstimate( RoP(ii),0.721,0.63,pc,ru(ii) );
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(ra(ii),ru(ii),RoP(ii),vacup,0.63,omega(ii),alpha(ii),gamma(ii,jj),pc );
    end
    %Fix omega at 0.93
    jj=4;
    parfor ii=1:ssize
        gamma(ii,jj) = gammaEstimate( RoP(ii),0.721,eps(ii),pc,ru(ii) );
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(ra(ii),ru(ii),RoP(ii),vacup,eps(ii),0.94,alpha(ii),gamma(ii,jj),pc );
    end
    %Fix ru at 2*10^(-4)
    jj=5;
    parfor ii=1:ssize
        gamma(ii,jj) = gammaEstimate( RoP(ii),0.721,eps(ii),pc,2*10^(-4) );
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(2*10^(-4).*fra(ii),2*10^(-4),RoP(ii),vacup,eps(ii),omega(ii),alpha(ii),gamma(ii,jj),pc );
    end
    %Fix fra at 5
    jj=6;
    parfor ii=1:ssize
        gamma(ii,jj) = gammaEstimate( RoP(ii),0.721,eps(ii),pc,ru(ii) );
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(ru(ii).*5,ru(ii),RoP(ii),vacup,eps(ii),omega(ii),alpha(ii),gamma(ii,jj),pc );
    end
    %Fix alpha
    jj=7;
    parfor ii=1:ssize
        gamma(ii,jj) = gammaEstimate( RoP(ii),0.721,eps(ii),pc,ru(ii) );
        [rho(ii,jj),minavgkappa(ii,jj),kappam(ii,jj)] = TableValues(ra(ii),ru(ii),RoP(ii),vacup,eps(ii),omega(ii),0.548,gamma(ii,jj),pc );
    end

    for jj=1:7
        f1=fopen(['Sensitivity to infection-pc=' num2str(pc) '-Case=' num2str(jj) '.txt'],'w');
        f2=fopen(['Density of prosocials-pc=' num2str(pc) '-Case=' num2str(jj) '.txt'],'w');
        f3=fopen(['MinKappa-pc=' num2str(pc) '-Case=' num2str(jj) '.txt'],'w');
        f4=fopen(['Kappa-pc=' num2str(pc) '-Case=' num2str(jj) '.txt'],'w');
        for ii=1:ssize
            if(~isinf(-min(kappam(ii,:))))
                fprintf(f1,'%32.31f \n',gamma(ii,jj));
                fprintf(f2,'%32.31f \n',rho(ii,jj));
                fprintf(f3,'%32.31f \n',max([min([minavgkappa(ii,jj),1]),0]));
                fprintf(f4,'%32.31f \n',max([min([kappam(ii,jj),1]),0]));
            end
        end
        fclose('all');
    end
 end
%delete(pobj);
% 
% for ii=2:7
%     subplot('Position',[xs(1) ys(pc) wt ht]);plot(linspace(min(min(gamma))*0.95,max(max(gamma))*1.05,100),ng(:,ii),'color',colors(ii,:),'LineWidth',2);hold on;
%     subplot('Position',[xs(2) ys(pc) wt ht]);plot(linspace(min(min(rho))*0.95,max(max(rho))*1.05,100),nr(:,ii),'color',colors(ii,:),'LineWidth',2);hold on;
%     subplot('Position',[xs(3) ys(pc) wt ht]);plot(linspace(min(min(kappam))*0.95,max(max(kappam))*1.05,100),nk(:,ii),'color',colors(ii,:),'LineWidth',2);hold on;
% end
% ii=1;
% subplot('Position',[xs(1) ys(pc) wt ht]);plot(linspace(min(min(gamma))*0.95,max(max(gamma))*1.05,100),ng(:,ii),'color',colors(ii,:),'LineWidth',2);hold off;
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',20);
% if(pc==3)
%     xlabel('Sensitivity to infection (\gamma)','Fontsize',24,'FontWeight','bold');
% end
% ylabel('Probability','Fontsize',24,'FontWeight','bold');
% subplot('Position',[xs(2) ys(pc) wt ht]);plot(linspace(min(min(rho))*0.95,max(max(rho))*1.05,100),nr(:,ii),'color',colors(ii,:),'LineWidth',2);hold off;
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',20);
% if(pc==3)
%     xlabel('Degree of prosociality (\rho)','Fontsize',24,'FontWeight','bold');
% end
% subplot('Position',[xs(3) ys(pc) wt ht]);plot(linspace(min(min(kappam))*0.95,max(max(kappam))*1.05,100),nk(:,ii),'color',colors(ii,:),'LineWidth',2);hold off;
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',20);
% if(pc==3)
%     xlabel('Quality of prosocial behavior (\kappa)','Fontsize',24,'FontWeight','bold');
% end
% V(:,1,pc)=1-var(gamma(:,2:7))./var(gamma(:,1));
% V(:,2,pc)=1-var(rho(:,2:7))./var(rho(:,1));
% V(:,3,pc)=1-var(kappam(:,2:7))./var(kappam(:,1));
% end
% 
% cn=[];
% for ii=2:7
% cn=[cn;colors(ii,:);colors(ii,:).*0.99;colors(ii,:).*0.98];
% end
% cn=[cn;[0 0 0];[0.01 0.01 0.01]; [0.02 0.02 0.02]];
% fHand = figure(10);
% aHand = axes('parent', fHand);
% hold(aHand, 'on')
% bar(-0.6, 0, 'parent', aHand, 'facecolor', [0 0 0],'barwidth',0.375);
% bar(-1, 0, 'parent', aHand, 'facecolor', [0.01 0.01 0.01],'barwidth',0.375);
% bar(-1.4, 0, 'parent', aHand, 'facecolor', [0.02 0.02 0.02],'barwidth',0.375);
% for ii=1:6
% bar(0.6+(ii-1)*1.5, V(ii,1,1), 'parent', aHand, 'facecolor', colors(ii+1,:),'barwidth',0.375);
% bar(1+(ii-1)*1.5, V(ii,1,2), 'parent', aHand, 'facecolor', colors(ii+1,:).*0.99,'barwidth',0.375);
% bar(1.4+(ii-1)*1.5, V(ii,1,3), 'parent', aHand, 'facecolor', colors(ii+1,:).*0.98,'barwidth',0.375);
% end
% ylabel('Relative decrease in variance','FontSize',28,'FontWeight','bold');
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',20,'XTick',[1:1.5:8.5],'XTicklabel','','YTick',[-1:0.1:1]);ylim([-1 1.02]);xlim([0 10.5]);
% set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
% [hx,hy] = format_ticks(gca,{'$R_0^P$','$\varepsilon$','$\omega$','$\psi$','$\xi$','$\alpha$'},[],[1:1.5:8.5],[],[],[],0.035,[]); set(hx,'FontSize',20);
% legend('Hill function','Max. prevelance','Exponential decay');
% legend('boxoff');
% applyhatch_pluscolor(gcf,'wxc',1,[1 0 1],cn)
% fHand = figure(12);
% aHand = axes('parent', fHand);
% hold(aHand, 'on')
% for ii=1:6
% bar(0.6+(ii-1)*1.5, V(ii,2,1), 'parent', aHand, 'facecolor', colors(ii+1,:),'barwidth',0.375);
% bar(1+(ii-1)*1.5, V(ii,2,2), 'parent', aHand, 'facecolor', colors(ii+1,:).*0.99,'barwidth',0.375);
% bar(1.4+(ii-1)*1.5, V(ii,2,3), 'parent', aHand, 'facecolor', colors(ii+1,:).*0.98,'barwidth',0.375);
% end
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',20,'XTick',[1:1.5:8.5],'XTicklabel','','YTick',[-1:0.1:1]);ylim([-1 1.02])
% set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
% [hx,hy] = format_ticks(gca,{'$R_0^P$','$\varepsilon$','$\omega$','$\psi$','$\xi$','$\alpha$'},[],[1:1.5:8.5],[],[],[],0.035,[]); set(hx,'FontSize',20);
% applyhatch_pluscolor(gcf,'wxc',1,[1 0 1],cn)
% fHand = figure(14);
% aHand = axes('parent', fHand);
% hold(aHand, 'on')
% for ii=1:6
% bar(0.6+(ii-1)*1.5, V(ii,3,1), 'parent', aHand, 'facecolor', colors(ii+1,:),'barwidth',0.375);
% bar(1+(ii-1)*1.5, V(ii,3,2), 'parent', aHand, 'facecolor', colors(ii+1,:).*0.99,'barwidth',0.375);
% bar(1.4+(ii-1)*1.5, V(ii,3,3), 'parent', aHand, 'facecolor', colors(ii+1,:).*0.98,'barwidth',0.375);
% end
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',20,'XTick',[1:1.5:8.5],'XTicklabel','','YTick',[-1:0.1:1]);ylim([-1 1.02])
% set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
% [hx,hy] = format_ticks(gca,{'$R_0^P$','$\varepsilon$','$\omega$','$\psi$','$\xi$','$\alpha$'},[],[1:1.5:8.5],[],[],[],0.035,[]); set(hx,'FontSize',20);
% applyhatch_pluscolor(gcf,'wxc',1,[1 0 1],cn)
