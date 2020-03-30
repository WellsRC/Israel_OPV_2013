Baselineparameters;
[rap rup]= Estrisk;

alpha=0.548;
RoP = RoPEstimate(rup,0.88,eps); %Set to the vaccine coverage of arab population 
pc=1;
[gamma] = gammaEstimate( RoP,0.721,eps,pc,rup );
[gamma2] = gammaEstimate( RoP,0.636,eps,pc,rup );
[rho minavgkappa kappa] = TableValues( rap,rup,RoP,0.721,eps,omega,alpha,gamma,pc);
[NEB1 NEB2 NEB3]=SolveNE(rap,rup,RoP,eps,kappa,omega,rho,gamma,alpha,pc);
VacupB=omega*alpha*rho*NEB1+omega*alpha*(1-rho)*NEB2+omega*(1-alpha)*NEB3;

[NEB1 NEB2 NEB3]=SolveNE(rap,rup,RoP,eps,kappa,omega,rho,gamma2,alpha,pc);
VacupB2=omega*alpha*rho*NEB1+omega*alpha*(1-rho)*NEB2+omega*(1-alpha)*NEB3;

alphav=linspace(0,1,101);
rhov=linspace(0,1,100);
VacupM=zeros(length(rhov),length(alphav));
VacupMK1=zeros(length(rhov),length(alphav));
VacupMR=zeros(length(rhov),length(alphav));
VacupMG=zeros(length(rhov),length(alphav));
%pobj=parpool(6);
for ii=1:length(rhov)
    parfor jj=1:length(alphav)
        [NE1 NE2 NE3]=SolveNE(rap,rup,RoP,eps,kappa,omega,rhov(ii),gamma,alphav(jj),pc);
        VacupM(ii,jj)=omega*alphav(jj)*rhov(ii)*NE1+omega*alphav(jj)*(1-rhov(ii))*NE2+omega*(1-alphav(jj))*NE3-VacupB;
        [NE1 NE2 NE3]=SolveNE(rap,rup,RoP,eps,minavgkappa,omega,rhov(ii),gamma,alphav(jj),pc);
        VacupMK1(ii,jj)=omega*alphav(jj)*rhov(ii)*NE1+omega*alphav(jj)*(1-rhov(ii))*NE2+omega*(1-alphav(jj))*NE3-VacupB;
        [NE1 NE2 NE3]=SolveNE(rup*1.01,rup,RoP,eps,kappa,omega,rhov(ii),gamma,alphav(jj),pc);
        VacupMR(ii,jj)=omega*alphav(jj)*rhov(ii)*NE1+omega*alphav(jj)*(1-rhov(ii))*NE2+omega*(1-alphav(jj))*NE3-VacupB;
        [NE1 NE2 NE3]=SolveNE(rap,rup,RoP,eps,kappa,omega,rhov(ii),gamma2,alphav(jj),pc);
        VacupG(ii,jj)=omega*alphav(jj)*rhov(ii)*NE1+omega*alphav(jj)*(1-rhov(ii))*NE2+omega*(1-alphav(jj))*NE3-VacupB2;
    end
end
% delete(pobj);
close all;
xs=[0.065 0.565];
ys=[0.575 0.08];
hei=0.33;
wid=0.37;
figure('units','normalized','outerposition',[0 0 1 1]);
minv=0.6;
x1=linspace(minv-VacupB,0,25);
x2=linspace(0,0.85-VacupB,25);
cv4=[x1(1:end-1) 0 x2(2:end)];

rmap=[103 0 13; 165 15 21; 203 24 29; 239 59 44; 251 106 74; 252 146 114; 252 187 161; 254 224 210; 255 245 240; 256 256 256]./256;
bmap=[256 256 256; 247 251 255; 198 219 239; 158 202 225; 107 174 214; 66 146 198; 33 113 181; 8 81 156; 8 48 107]./256;
ncol=64;
f=find(linspace(min(cv4),max(cv4),ncol)>0);
g=find(linspace(min(cv4),max(cv4),ncol)<=0);
test=linspace(min(cv4),max(cv4),ncol);
if(abs(test(g(end)))<abs(test(f(1))))
h=g(end);
else
h=f(1);
end
if(h<2)
    idx3 = linspace(0,1,size(bmap,1));
    idx4 = linspace(0,1,64);
    map = [interp1(idx3,bmap,idx4)] ;
else
    idx1 = linspace(0,1,size(rmap,1));
    idx2 = linspace(0,1,h);    
    idx3 = linspace(0,1,size(bmap,1));
    idx4 = linspace(0,1,64-h);
    map = [interp1(idx1,rmap,idx2);interp1(idx3,bmap,idx4)] ;
end
colormap(map)
subplot('Position',[xs(1) ys(2) wid hei])
contourf(alphav,rhov,VacupMK1,cv4,'LineStyle','none');caxis([minv-VacupB 0.85-VacupB]);
box off;
set(gca,'LineWidth',2,'TickDir','out','Xlim',[0 1],'Ylim',[0 1],'XMinorTick','on','YMinorTick','on','XTick',[0:0.1:1],'YTick',[0:0.1:1],'Fontsize',18);
xlabel('Level of comprehension (\alpha)','Fontsize',22,'Fontweight','bold');
ylabel({'Prevalence of','prosociality (\rho)'},'Fontsize',22,'Fontweight','bold');
title({['Low strength of prosocial'],['behavior (\kappa=' num2str(round(minavgkappa,2)) ')']},'Fontsize',22,'Fontweight','bold');
h=colorbar('EastOutside');
set(h,'LineWidth',2,'YTick',[minv-VacupB:0.05:0.85-VacupB],'YTicklabel',{'60%','65%','70%','75%','80%','85%'},'Fontsize',18,'Position',[xs(1)+wid+0.0025,ys(2),0.011204481792717,0.330451488952931])
ylabel(h,{'Vaccination coverage (\nu)'},'Fontsize',22,'Fontweight','bold','rotation',270,'position',[5.355555443536657,-0.0]);
h=text(-0.163465909090909,1.16,'C');
set(h,'Fontsize',34,'Fontweight','bold');

subplot('Position',[xs(2) ys(1) wid hei])
contourf(alphav,rhov,VacupMR,cv4,'LineStyle','none');caxis([minv-VacupB 0.85-VacupB]);
box off;
set(gca,'LineWidth',2,'TickDir','out','Xlim',[0 1],'Ylim',[0 1],'XMinorTick','on','YMinorTick','on','XTick',[0:0.1:1],'YTick',[0:0.1:1],'Fontsize',18);
xlabel('Level of comprehension (\alpha)','Fontsize',22,'Fontweight','bold');
ylabel({'Prevalence of','prosociality (\rho)'},'Fontsize',22,'Fontweight','bold');
title({'Relative risk comparable','between groups (r_A\approx r_U)'},'Fontsize',22,'Fontweight','bold');
h=text(-0.163465909090909,1.16,'B');
set(h,'Fontsize',34,'Fontweight','bold');
h=colorbar('EastOutside');
set(h,'LineWidth',2,'YTick',[minv-VacupB:0.05:0.85-VacupB],'YTicklabel',{'60%','65%','70%','75%','80%','85%'},'Fontsize',18,'Position',[xs(2)+wid+0.0025,ys(1),0.011204481792717,0.330451488952931])

ylabel(h,{'Vaccination coverage (\nu)'},'Fontsize',22,'Fontweight','bold','rotation',270,'position',[5.355555443536657,-0.0]);

subplot('Position',[xs(1) ys(1) wid hei])
contourf(alphav,rhov,VacupM,cv4,'LineStyle','none');caxis([minv-VacupB 0.85-VacupB]);
box off;
set(gca,'LineWidth',2,'TickDir','out','Xlim',[0 1],'Ylim',[0 1],'XMinorTick','on','YMinorTick','on','XTick',[0:0.1:1],'YTick',[0:0.1:1],'Fontsize',18);
xlabel('Level of comprehension (\alpha)','Fontsize',22,'Fontweight','bold');
ylabel({'Prevalence of','prosociality (\rho)'},'Fontsize',22,'Fontweight','bold');
title({'Baseline',['(\kappa=' num2str(round(kappa,2)) ')']},'Fontsize',22,'Fontweight','bold');
h=text(-0.163465909090909,1.16,'A');
set(h,'Fontsize',34,'Fontweight','bold');
h=colorbar('EastOutside');
set(h,'LineWidth',2,'YTick',[minv-VacupB:0.05:0.85-VacupB],'YTicklabel',{'60%','65%','70%','75%','80%','85%'},'Fontsize',18,'Position',[xs(1)+wid+0.0025,ys(1),0.011204481792717,0.330451488952931])

ylabel(h,{'Vaccination coverage (\nu)'},'Fontsize',22,'Fontweight','bold','rotation',270,'position',[5.355555443536657,-0.0]);



figure('units','normalized','outerposition',[0 0 1 1]);
minv=0.5;
x1=linspace(minv-VacupB2,0,25);
x2=linspace(0,0.7-VacupB2,25);
cv4=[x1(1:end-1) 0 x2(2:end)];

rmap=[103 0 13; 165 15 21; 203 24 29; 239 59 44; 251 106 74; 252 146 114; 252 187 161; 254 224 210; 255 245 240; 256 256 256]./256;
bmap=[256 256 256; 247 251 255; 198 219 239; 158 202 225; 107 174 214; 66 146 198; 33 113 181; 8 81 156; 8 48 107]./256;
ncol=64;
f=find(linspace(min(cv4),max(cv4),ncol)>0);
g=find(linspace(min(cv4),max(cv4),ncol)<=0);
test=linspace(min(cv4),max(cv4),ncol);
if(abs(test(g(end)))<abs(test(f(1))))
h=g(end);
else
h=f(1);
end
if(h<2)
    idx3 = linspace(0,1,size(bmap,1));
    idx4 = linspace(0,1,64);
    map = [interp1(idx3,bmap,idx4)] ;
else
    idx1 = linspace(0,1,size(rmap,1));
    idx2 = linspace(0,1,h);    
    idx3 = linspace(0,1,size(bmap,1));
    idx4 = linspace(0,1,64-h);
    map = [interp1(idx1,rmap,idx2);interp1(idx3,bmap,idx4)] ;
end
colormap(map)

subplot('Position',[xs(2) ys(2) wid hei])
contourf(alphav,rhov,VacupG,cv4,'LineStyle','none');caxis([minv-VacupB2 0.7-VacupB2]);
box off;
set(gca,'LineWidth',2,'TickDir','out','Xlim',[0 1],'Ylim',[0 1],'XMinorTick','on','YMinorTick','on','XTick',[0:0.1:1],'YTick',[0:0.1:1],'Fontsize',18);
xlabel('Level of comprehension (\alpha)','Fontsize',22,'Fontweight','bold');
ylabel({'Prevalence of','prosociality (\rho)'},'Fontsize',22,'Fontweight','bold');
title({['Lower perceived likelihood of infection'],['(\gamma=' num2str(round(gamma2,3)) ')']},'Fontsize',22,'Fontweight','bold');
h=text(-0.163465909090909,1.16,'D');
set(h,'Fontsize',34,'Fontweight','bold');
h=colorbar('EastOutside');
set(h,'LineWidth',2,'YTick',[minv-VacupB2:0.05:0.7-VacupB2],'YTicklabel',{'50%','55%','60%','65%','70%','75%'},'Fontsize',18,'Position',[xs(2)+wid+0.0025,ys(2),0.011204481792717,0.330451488952931])

ylabel(h,{'Vaccination coverage (\nu)'},'Fontsize',22,'Fontweight','bold','rotation',270,'position',[5.355555443536657,-0.039999904632569]);

%clear;