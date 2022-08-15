% This file is used to plot the figure 5
% Author: Lei Zhang
% Last modified: 2022-06-27

clearvars;
close all;
clc;

%%
% read the parameter of 23 provinces
tableSD = readtable('province.csv');
location = tableSD.location;

alpha = tableSD.alpha;
h1 = tableSD.h0;
th = tableSD.t_h;
h0 = tableSD.H1;
duration_cured = tableSD.t_cure;

%% fit the alpha and th

% Remove provinces that do not have a full evolution cycle
province_remove = {'Anhui','Beijing','Guangxi','Guizhou','Shandong','Shanghai'};
indx = ismember(location,province_remove);
location_removed = location(~indx);
alpha_removed = alpha(~indx);
h1_removed = h1(~indx);
h0_removed = h0(~indx);
th_removed = th(~indx);
duration_cured_removed = duration_cured(~indx);

% Distinguish two patterns
pattern_1 = {'Henan','Yunnan','Hebei','Fujian','Jiangxi','Jiangsu','Shaanxi','Tianjin'};
pattern_2 = {'Hainan','Hunan','Chongqing','Shanxi','Zhejiang','Heilongjiang',...
    'Guangdong','Sichuan','Liaoning'};

indx1 = ismember(location_removed,pattern_1);
location_removed1 = location_removed(indx1);
alpha_removed1 = alpha_removed(indx1);
h1_removed1 = h1_removed(indx1);
h0_removed1 = h0_removed(indx1);
th_removed1 = th_removed(indx1);
duration_cured_removed1 = duration_cured_removed(indx1);

h1_h0_1 = mean(h1_removed1./h0_removed1);

mean_duration1 = mean(duration_cured_removed1);
duration1_25   = quantile(duration_cured_removed1,0.25);
duration1_75   = quantile(duration_cured_removed1,0.75);

disp(['The duration of pattern 1 is ',num2str(mean_duration1),...
    '(IQR, ',num2str(duration1_25),'，',num2str(duration1_75),')'])

indx2 = ismember(location_removed,pattern_2);
location_removed2 = location_removed(indx2);
alpha_removed2 = alpha_removed(indx2);
h1_removed2 = h1_removed(indx2);
h0_removed2 = h0_removed(indx2);
th_removed2 = th_removed(indx2);
duration_cured_removed2 = duration_cured_removed(indx2);

h1_h0_2 = mean(h1_removed2./h0_removed2);

ratio_h1_h0 = h1_h0_1 / h1_h0_2;

mean_duration2 = mean(duration_cured_removed2);
duration2_25   = quantile(duration_cured_removed2,0.25);
duration2_75   = quantile(duration_cured_removed2,0.75);

disp(['The duration of pattern 2 is ',num2str(mean_duration2),...
    '(IQR, ',num2str(duration2_25),'，',num2str(duration2_75),')'])

disp(['The ratio of h1_h0 of pattern 1 and 2  is ',num2str(ratio_h1_h0)])

% fit the alpha-th of pattern1
[fitresult_alpha_th1,gof_alpha_th1] = fit_alpha_th(alpha_removed1,th_removed1);
x1_alpha_th1 = fitresult_alpha_th1.x1;
x2_alpha_th1 = fitresult_alpha_th1.x2;

alpha_pre = 0.05:0.001:0.25;
th_pre1 = x1_alpha_th1 .* alpha_pre .^ x2_alpha_th1;

r2_alpha_th1 = gof_alpha_th1.rsquare;


% fit the alpha-th of pattern2
[fitresult_alpha_th2,gof_alpha_th2] = fit_alpha_th(alpha_removed2,th_removed2);
x1_alpha_th2 = fitresult_alpha_th2.x1;
x2_alpha_th2 = fitresult_alpha_th2.x2;

th_pre2 = x1_alpha_th2 .* alpha_pre .^ x2_alpha_th2;

r2_alpha_th2 = gof_alpha_th2.rsquare;


%% fit the h0 and h1

% fit the h0-h1 of pattern1
[fitresult_h0_h1_1,gof_h0_h1_1] = fit_h0_h1(h0_removed1,h1_removed1);
x1_h0_h1_1 = fitresult_h0_h1_1.x1;
x2_h0_h1_1 = fitresult_h0_h1_1.x2;

h0_pre = 0.0001:0.001:1;
h1_pre_1 = x1_h0_h1_1 .* h0_pre + x2_h0_h1_1;

r2_h0_h1_1 = gof_h0_h1_1.rsquare;

% fit the h0-h1 of pattern2
[fitresult_h0_h1_2,gof_h0_h1_2] = fit_h0_h1(h0_removed2,h1_removed2);
x1_h0_h1_2 = fitresult_h0_h1_2.x1;
x2_h0_h1_2 = fitresult_h0_h1_2.x2;

h1_pre_2 = x1_h0_h1_2 .* h0_pre + x2_h0_h1_2;

r2_h0_h1_2 = gof_h0_h1_2.rsquare;




%% Plot the results

gcf = figure('position',[382,107,1186,809]);
t = tiledlayout(21,2,'TileSpacing','tight','Padding','compact');

gca = nexttile([9,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on

plot(alpha_removed1,th_removed1,'o',...
    'color',[0.80392	0.33333	0.33333],...
    'LineWidth',1.5,'markersize',8,'markerfacecolor',[0.80392	0.33333	0.33333])
plot(alpha_removed2,th_removed2,'o',...
    'color',[0.27451	0.5098	0.70588],...
    'LineWidth',1.5,'markersize',8,'markerfacecolor',[0.27451	0.5098	0.70588])

plot(alpha_pre,th_pre1,'-','color',[0.80392	0.33333	0.33333],'LineWidth',2)
plot(alpha_pre,th_pre2,'k-','color',[0.27451	0.5098	0.70588],'LineWidth',2)

xlabel('growth rate of healing degree ($\alpha$)',...
    'FontName','Times New Roman','FontSize',16,'interpreter','latex')
ylabel('inflection point ($t_h$)',...
    'FontName','Times New Roman','FontSize',16,'interpreter','latex')

legend({'pattern 1','pattern 2',...
    ['y = ',sprintf('%0.2f*x^{%0.2f}',x1_alpha_th1,x2_alpha_th1)],...
    ['y = ',sprintf('%0.2f*x^{%0.2f}',x1_alpha_th2,x2_alpha_th2)]},...
    'location','northeast','FontSize',16,'numcolumns',1)
legend('boxoff')

text(0.06,27,['R^2 = ' ' ' num2str(roundn(r2_alpha_th1,-2))],...
    'color',[0.80392	0.33333	0.33333],...
    'LineWidth',1,'FontSize',14)
text(0.06,23.6,['R^2 = ' ' ' num2str(roundn(r2_alpha_th2,-2))],...
    'color',[0.27451	0.5098	0.70588],...
    'LineWidth',1,'FontSize',14)

title(['$(a) \quad \alpha$ ' '---' ' $t_h$'],'Position',[0.15,5,0],...
    'FontSize',18,'interpreter','latex')

set(gca,'FontSize',16)
box on
xlim([0.05,0.25])
ylim([20,60])


%%
gca = nexttile([9,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on

plot(h0_removed1,h1_removed1,'o',...
    'color',[0.80392	0.33333	0.33333],...
    'LineWidth',1.5,'markersize',8,'markerfacecolor',[0.80392	0.33333	0.33333])
plot(h0_removed2,h1_removed2,'o',...
    'color',[0.27451	0.5098	0.70588],...
    'LineWidth',1.5,'markersize',8,'markerfacecolor',[0.27451	0.5098	0.70588])

plot(h0_pre,h1_pre_1,'-','color',[0.80392	0.33333	0.33333],'LineWidth',2)
plot(h0_pre,h1_pre_2,'k-','color',[0.27451	0.5098	0.70588],'LineWidth',2)

xlabel('initial healing degree ($H_0$)',...
    'FontName','Times New Roman','FontSize',16,'interpreter','latex')
ylabel('saturation healing degree ($H_1$)',...
    'FontName','Times New Roman','FontSize',16,'interpreter','latex')

legend({'pattern 1','pattern 2',...
    ['y = ',sprintf('%0.2f*x+%0.2f',x1_h0_h1_1,x2_h0_h1_1)],...
    ['y = ',sprintf('%0.2f*x+%0.2f',x1_h0_h1_2,x2_h0_h1_2)]},...
    'location','northwest','FontSize',16,'numcolumns',1)
legend('boxoff')

text(0.0068,0.8,['R^2 = ' ' ' num2str(roundn(r2_h0_h1_1,-2))],...
    'color',[0.80392	0.33333	0.33333],...
    'LineWidth',1,'FontSize',14)
text(0.0068,0.69,['R^2 = ' ' ' num2str(roundn(r2_h0_h1_2,-2))],...
    'color',[0.27451	0.5098	0.70588],...
    'LineWidth',1,'FontSize',14)

title(['$(b) \quad H_0$' '---' '$H_1$'],'Position',[0.0045,-0.45,0],...
    'FontSize',18,'interpreter','latex')

set(gca,'FontSize',16)
box on
xlim([0,0.009])
ylim([0,1.2])
%%

gca = nexttile([1,2]);
gca.Visible = 'off';


%%
gca = nexttile([9,2]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on

bar_1 = length(duration_cured_removed1);
bar_2 = length(duration_cured_removed2);

b1 = bar(1:1:bar_1,sort(duration_cured_removed1),0.5,...
    'facecolor',[0.80392	0.33333	0.33333],...
    'edgecolor',[0.80392	0.33333	0.33333]);
b2 = bar(bar_1+1:1:bar_1+bar_2,sort(duration_cured_removed2),0.5,...
    'facecolor',[0.27451	0.5098	0.70588],...
    'edgecolor',[0.27451	0.5098	0.70588]);

xtips1 = b1.XEndPoints;
ytips1 = b1.YEndPoints;
labels1 = string(b1.YData);

xtips2 = b2.XEndPoints;
ytips2 = b2.YEndPoints;
labels2 = string(b2.YData);

text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',14)
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',14)

ylabel('duration (days)',...
    'FontName','Times New Roman','FontSize',16,'interpreter','latex')

xticks(1:1:bar_1+bar_2)
xticklabels([pattern_1,pattern_2]);

xlim([0,bar_1+bar_2+1])
ylim([0,60])

set(gca,'fontsize',16)
box on
set(gca,'Layer','top');

title(['$(a) \quad \alpha$ ' '---' ' $t_h$'],'Position',[0.15,5,0],...
    'FontSize',18,'interpreter','latex')


title('(c)\quad duration of cured cases','Position',[9,-27,0],...
    'FontName','Times New Roman','FontSize',18,'interpreter','latex')

%%

gca = nexttile([2,2]);
gca.Visible = 'off';
%%

function [fitresult,gof] = fit_alpha_th(alpha,th)

[xData,yData] = prepareCurveData(alpha,th);

% Set up fittype and options.
ft = fittype( 'x1.*x.^x2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [7.06265602545247 -0.75507766185979];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end


function [fitresult,gof] = fit_h0_h1(alpha,th)

[xData,yData] = prepareCurveData(alpha,th);

% Set up fittype and options.
ft = fittype( 'x1*x+x2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.421761282626275 0.915735525189067];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end
