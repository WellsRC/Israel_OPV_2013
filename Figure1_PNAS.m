close all;
xs=[0.0638 0.3788 0.705];
ys=[0.6 0.1];
hei=0.30;
wid=0.28;
alpha=0.548;
rho=0.655/0.94;
omega=0.94;
n=1;

figure('units','normalized','outerposition',[0 0 1 1]);
p=load('721VacCoverage-HeterogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt');
p2=load('721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt');
sizeg=100;
gamma=[linspace(0.05,0.2,sizeg);];
xt=zeros(size(gamma));
yb=prctile(p',97.5);
yt=zeros(size(yb));
for ii=1:sizeg
    xt(ii)=gamma(end-(ii-1));
    yt(ii)=yb(end-(ii-1));
end
x1=[gamma xt];
ya=prctile(p',2.5);
y1=[ya yt];
subplot('Position',[xs(1) ys(1) wid hei]);p3=patch(x1,y1,[118 42 131]./256,'FaceAlpha',.25);
    set(p3,'EdgeColor',[118 42 131]./256,'LineStyle','none'); hold on
    plot(gamma,median(p'),'color',[118 42 131]./256,'LineWidth',2);
    yb=prctile(p2',97.5);
yt=zeros(size(yb));
for ii=1:sizeg
    yt(ii)=yb(end-(ii-1));
end
ya=prctile(p2',2.5);
y1=[ya yt];
p3=patch(x1,y1,[1 102 94]./256,'FaceAlpha',.25);
    set(p3,'EdgeColor',[1 102 94]./256,'LineStyle','none'); 
    plot(gamma,median(p2'),'color',[1 102 94]./256,'LineWidth',2);hold off
title('National vaccination coverage (\nu)','FontSize',22,'FontWeight','bold');
box off;
set(gca,'Linewidth',2,'tickdir','out','Fontsize',20,'XTick',[0.05:0.025:0.2],'YTick',[0:0.1:1],'Xminortick','on','Yminortick','on');ylim([0 1]);xlim([0.05 0.2]);
ylabel({'Vaccination','coverage (\nu)'},'FontSize',22,'FontWeight','bold');
xlabel({'Sensitivity to infection (\gamma)'},'FontSize',22,'FontWeight','bold');
text(0.016730337078652,1.070945945945946,'A','Fontsize',34,'FontWeight','bold');

p=omega.*(rho.*load('ProStrat-721VacCoverage-HeterogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt')+(1-rho).*load('IndStrat-721VacCoverage-HeterogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt'));
p2=omega.*(rho.*load('ProStrat-721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt')+(1-rho).*load('IndStrat-721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt'));

sizeg=100;
gamma=[linspace(0.05,0.2,sizeg);];
xt=zeros(size(gamma));
yb=prctile(p',97.5);
yt=zeros(size(yb));
for ii=1:sizeg
    xt(ii)=gamma(end-(ii-1));
    yt(ii)=yb(end-(ii-1));
end
x1=[gamma xt];
ya=prctile(p',2.5);
y1=[ya yt];
subplot('Position',[xs(2) ys(1) wid hei]);p3=patch(x1,y1,[118 42 131]./256,'FaceAlpha',.25);
    set(p3,'EdgeColor',[118 42 131]./256,'LineStyle','none'); hold on
    plot(gamma,median(p'),'color',[118 42 131]./256,'LineWidth',2);
    yb=prctile(p2',97.5);
yt=zeros(size(yb));
for ii=1:sizeg
    yt(ii)=yb(end-(ii-1));
end
ya=prctile(p2',2.5);
y1=[ya yt];
p3=patch(x1,y1,[1 102 94]./256,'FaceAlpha',.25);
    set(p3,'EdgeColor',[1 102 94]./256,'LineStyle','none'); 
    plot(gamma,median(p2'),'color',[1 102 94]./256,'LineWidth',2);hold off
title('Vaccination coverage among aware','FontSize',22,'FontWeight','bold');
box off;
set(gca,'Linewidth',2,'tickdir','out','Fontsize',20,'XTick',[0.05:0.025:0.2],'YTick',[0:0.1:1],'Xminortick','on','Yminortick','on');ylim([0 1]);xlim([0.05 0.2]);
xlabel({'Sensitivity to infection (\gamma)'},'FontSize',22,'FontWeight','bold');
text(0.035,1.070945945945946,'B','Fontsize',34,'FontWeight','bold');

p=omega.*(load('UnawareStrat-721VacCoverage-HeterogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt'));
p2=omega.*(load('UnawareStrat-721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt'));

sizeg=100;
gamma=[linspace(0.05,0.2,sizeg);];
xt=zeros(size(gamma));
yb=prctile(p',97.5);
yt=zeros(size(yb));
for ii=1:sizeg
    xt(ii)=gamma(end-(ii-1));
    yt(ii)=yb(end-(ii-1));
end
x1=[gamma xt];
ya=prctile(p',2.5);
y1=[ya yt];
subplot('Position',[xs(3) ys(1) wid hei]);p3=patch(x1,y1,[118 42 131]./256,'FaceAlpha',.25);
    set(p3,'EdgeColor',[118 42 131]./256,'LineStyle','none'); hold on
    plot(gamma,median(p'),'color',[118 42 131]./256,'LineWidth',2);
    yb=prctile(p2',97.5);
yt=zeros(size(yb));
for ii=1:sizeg
    yt(ii)=yb(end-(ii-1));
end
ya=prctile(p2',2.5);
y1=[ya yt];
p3=patch(x1,y1,[1 102 94]./256,'FaceAlpha',.25);
    set(p3,'EdgeColor',[1 102 94]./256,'LineStyle','none'); 
    plot(gamma,median(p2'),'color',[1 102 94]./256,'LineWidth',2);hold off
title('Vaccination coverage among unaware','FontSize',22,'FontWeight','bold');
box off;
set(gca,'Linewidth',2,'tickdir','out','Fontsize',20,'XTick',[0.05:0.025:0.2],'YTick',[0:0.1:1],'Xminortick','on','Yminortick','on');ylim([0 1]);xlim([0.05 0.2]);
xlabel({'Sensitivity to infection (\gamma)'},'FontSize',22,'FontWeight','bold');
text(0.035,1.070945945945946,'C','Fontsize',34,'FontWeight','bold');
text(0.142,0.535,'Heterogeneous','FontSize',22,'FontWeight','bold','color',[118 42 131]./256);
text(0.142,0.65,'Homogeneous','FontSize',22,'FontWeight','bold','color',[1 102 94]./256);

p=omega.*(load('ProStrat-721VacCoverage-HeterogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt'));
p2=omega.*(load('ProStrat-721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt'));
sizeg=100;
gamma=[linspace(0.05,0.2,sizeg);];
xt=zeros(size(gamma));
yb=prctile(p',97.5);
yt=zeros(size(yb));
for ii=1:sizeg
    xt(ii)=gamma(end-(ii-1));
    yt(ii)=yb(end-(ii-1));
end
x1=[gamma xt];
ya=prctile(p',2.5);
y1=[ya yt];
subplot('Position',[xs(1) ys(2) wid hei]);p3=patch(x1,y1,[118 42 131]./256,'FaceAlpha',.25);
    set(p3,'EdgeColor',[118 42 131]./256,'LineStyle','none'); hold on
    plot(gamma,median(p'),'color',[118 42 131]./256,'LineWidth',2);
    yb=prctile(p2',97.5);
yt=zeros(size(yb));
for ii=1:sizeg
    yt(ii)=yb(end-(ii-1));
end
ya=prctile(p2',2.5);
y1=[ya yt];
p3=patch(x1,y1,[1 102 94]./256,'FaceAlpha',.25);
    set(p3,'EdgeColor',[1 102 94]./256,'LineStyle','none'); 
    plot(gamma,median(p2'),'color',[1 102 94]./256,'LineWidth',2);hold off
title({'Vaccination coverage among','prosocials'},'FontSize',22,'FontWeight','bold');
box off;
set(gca,'Linewidth',2,'tickdir','out','Fontsize',20,'XTick',[0.05:0.025:0.2],'YTick',[0:0.1:1],'Xminortick','on','Yminortick','on');ylim([0 1]);xlim([0.05 0.2]);
ylabel({'Vaccination','coverage (\nu)'},'FontSize',22,'FontWeight','bold');
xlabel({'Sensitivity to infection (\gamma)'},'FontSize',22,'FontWeight','bold');
text(0.016730337078652,1.189,'D','Fontsize',34,'FontWeight','bold');


p=omega.*(load('IndStrat-721VacCoverage-HeterogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt'));
p2=omega.*(load('IndStrat-721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt'));
sizeg=100;
gamma=[linspace(0.05,0.2,sizeg);];
xt=zeros(size(gamma));
yb=prctile(p',97.5);
yt=zeros(size(yb));
for ii=1:sizeg
    xt(ii)=gamma(end-(ii-1));
    yt(ii)=yb(end-(ii-1));
end
x1=[gamma xt];
ya=prctile(p',2.5);
y1=[ya yt];
subplot('Position',[xs(2) ys(2) wid hei]);p3=patch(x1,y1,[118 42 131]./256,'FaceAlpha',.25);
    set(p3,'EdgeColor',[118 42 131]./256,'LineStyle','none'); hold on
    plot(gamma,median(p'),'color',[118 42 131]./256,'LineWidth',2);
    yb=prctile(p2',97.5);
yt=zeros(size(yb));
for ii=1:sizeg
    yt(ii)=yb(end-(ii-1));
end
ya=prctile(p2',2.5);
y1=[ya yt];
p3=patch(x1,y1,[1 102 94]./256,'FaceAlpha',.25);
    set(p3,'EdgeColor',[1 102 94]./256,'LineStyle','none'); 
    plot(gamma,median(p2'),'color',[1 102 94]./256,'LineWidth',2);hold off
title({'Vaccination coverage among','individualist'},'FontSize',22,'FontWeight','bold');
box off;
set(gca,'Linewidth',2,'tickdir','out','Fontsize',20,'XTick',[0.05:0.025:0.2],'YTick',[0:0.1:1],'Xminortick','on','Yminortick','on');ylim([0 1]);xlim([0.05 0.2]);
xlabel({'Sensitivity to infection (\gamma)'},'FontSize',22,'FontWeight','bold');
text(0.035,1.189,'E','Fontsize',34,'FontWeight','bold');

p=load('721VacCoverage-HeterogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt');
p2=load('721VacCoverage-HomogeneousRisk-MinKappa-ProsocialAware-rho=0.69681-pc=1-1.txt');
sigeg=100;
gamma=[linspace(0.05,0.2,sizeg);];
LH=zeros(sizeg,2);
for ii=1:sizeg
            temp=p(ii,:);
            n2=hist(temp,[0:0.01:1]);
            n2=n2./sum(n2);
            f=find([0:0.01:1]==0.72);
            LH(ii,1)=n2(f);
            temp=p2(ii,:);
            n2=hist(temp,[0:0.01:1]);
            n2=n2./sum(n2);
            f=find([0:0.01:1]==0.72);
            LH(ii,2)=n2(f);
end
A1=sum(LH(2:end,1).*(gamma(2)-gamma(1)));
A2=sum(LH(2:end,2).*(gamma(2)-gamma(1)));
subplot('Position',[xs(3) ys(2) wid hei]);
plot(gamma,LH(:,2),'color',[1 102 94]./256,'LineWidth',2); hold on;ylim([0,0.8]);
plot(gamma,LH(:,1),'color',[118 42 131]./256,'LineWidth',2); 
text(0.15006191369606,0.770675675675676,['AUC=' num2str(round(A2,4))],'Fontsize',22,'FontWeight','bold','color',[1 102 94]./256);
text(0.15006191369606,0.687,['AUC=' num2str(round(A1,4))],'Fontsize',22,'FontWeight','bold','color',[118 42 131]./256);
set(gca,'Linewidth',2,'tickdir','out','Fontsize',20,'XTick',[0.05:0.025:0.2],'YTick',[0:0.1:0.8],'Xminortick','on','Yminortick','on');xlim([0.05 0.2]);
xlabel({'Sensitivity to infection (\gamma)'},'FontSize',22,'FontWeight','bold');
text(0.024305816135084,0.952,'F','Fontsize',34,'FontWeight','bold');
box off;
ylabel('Likelihood','FontSize',22,'FontWeight','bold');
title({'Likelihood of 72%', 'national vaccination coverage'},'FontSize',22,'FontWeight','bold');
A1/A2
clear;