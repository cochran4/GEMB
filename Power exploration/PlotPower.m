function PlotPower
close all
% Effect sizes
delta = 2*10^(-4);

% Color levels
lvls1 = 0.5:0.1:3.5;
lvls2 = 0.5:1:4;


SubPlotPower(1000,10^(-16),10^(-16),lvls1)
title({'Rejection regions (area = 0.05)','under null hypothesis'},'FontWeight','normal')
ylabel('Gene 2 rank / n','FontSize',11)
xlabel('Gene 1 rank / n','FontSize',11)
legend({'Single','Gene set','Weighted gene set'})
legend boxoff
set(gcf,'Position',[50 50 400 300])
%print('fig_reject','-dpng','-r1000')
print('fig_reject','-dsvg','-r1000')

% Plot 3 figures
figure
subplot(2,2,1)
SubPlotPower(1000,delta,delta,lvls1)
title('k=1000, d_1 = d_2 = 0.0002','FontWeight','normal')
ylabel('Gene 2 rank / n','FontSize',11)

subplot(2,2,2)
SubPlotPower(1000,delta/2,3*delta/2,lvls1)
title('k=1000, d_1 = 0.0001, d_2 = 0.0003','FontWeight','normal')

subplot(2,2,3)
SubPlotPower(10000,delta,delta,lvls2)
title('k=10000, d_1 = d_2 = 0.0002','FontWeight','normal')
xlabel('Gene 1 rank / n','FontSize',11)
ylabel('Gene 2 rank / n','FontSize',11)

subplot(2,2,4)
SubPlotPower(10000,delta/2,3*delta/2,lvls2)
title('k=10000, d_1 = 0.0001, d_2 = 0.0003','FontWeight','normal')
xlabel('Gene 1 rank / n','FontSize',11)


set(gcf,'Position',[50 50 800 600])
print('fig_power','-dsvg','-r1000')


function SubPlotPower(k,effect1,effect2,lvls)


% Grid of pvalues
[p2,p1] = ndgrid(0.01:0.01:0.99);

% Density
nu = 10;
logdensity1 = log( ncfpdf(finv(1-p1,nu-1,k-nu),nu-1,k-nu,effect1*k) ) - log( fpdf(finv(1-p1,nu-1,k-nu),nu-1,k-nu) ); % (2*norminv(1-p1/2)-effect1*sqrt(n))*effect1*sqrt(n); 
logdensity2 = log( ncfpdf(finv(1-p2,nu-1,k-nu),nu-1,k-nu,effect2*k) ) - log( fpdf(finv(1-p2,nu-1,k-nu),nu-1,k-nu) ); %(2*norminv(1-p2/2)-effect2*sqrt(n))*effect2*sqrt(n); 

% Joint density
joint = exp(logdensity1 + logdensity2);


% Plot
Colors
PrettyFig
hold on

%patch(0.0253*[1; 1; 0],[0.0253 1],'Color',[1 1 1]*0.5,'LineWidth',1.25)
%patch([0.0253 1],0.0253*[1 1],'Color',[1 1 1]*0.5,'LineWidth',1.25)
pcrit = 0.0253;
patch([0 1 1 pcrit pcrit 0],[0 0 pcrit pcrit 1 1],clr(1,:),'FaceAlpha',0.2)
patch(sqrt(0.1)*[0 0 1],sqrt(0.1)*[0 1 0],clr(2,:),'FaceAlpha',0.2)
patch(3*sqrt(0.1/3)*[0 0 1],sqrt(0.1/3)*[0 1 0],clr(4,:),'FaceAlpha',0.2)



N = 10^4;
[C,h] = contour(p1,p2,joint,lvls,'LineWidth',1.25);
contourcmap('autumn',lvls)
clabel(C,h,'LabelSpacing',350,'Margin',4,'FontSize',9.5)
set(gca,'XLim',[0 1])
set(gca,'YLim',[0 1])
