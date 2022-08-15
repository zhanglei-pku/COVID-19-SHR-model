% This file is used to plot the figure 6
% Author: Lei Zhang
% Last modified: 2022-06-27

clear
close all
clc

inputarea = 'Beijing';
% read the parameters
tableBJ = readtable('Beijing.csv');
s1 = tableBJ.s1;
gamma = tableBJ.gamma;
alpha = tableBJ.alpha;
cum_cases = tableBJ.cum_cases;

% Call the fitting function
% the first wave
[DataTobeFitted1,data_pre1,N_1,time_pre1,...
    date_cum_confirmed1,date_cum_recovered1,...
    date_new_confirmed1,date_new_recovered1,...
    date_quarantined1] = func_fit_1_2wave(1,inputarea);

gap_cum_1 = datenum(date_cum_recovered1(1)-date_cum_confirmed1(1))+1;
gap_new_1 = datenum(date_new_recovered1(1)-date_cum_confirmed1(1))+1;

% the second wave
[DataTobeFitted2,data_pre2,N_2,time_pre2,...
    date_cum_confirmed2,date_cum_recovered2,...
    date_new_confirmed2,date_new_recovered2,...
    date_quarantined2] = func_fit_1_2wave(2,inputarea);

gap_cum_2 = datenum(date_cum_recovered2(1)-date_cum_confirmed2(1))+1;
gap_new_2 = datenum(date_new_recovered2(1)-date_cum_confirmed2(1))+1;

% 
Data_all = {DataTobeFitted1,DataTobeFitted2};
NN = {N_1,N_2};
data_all = {data_pre1,data_pre2};
time_all = {time_pre1,time_pre2};
date_cum_confirmed_all = {date_cum_confirmed1,date_cum_confirmed2};
date_cum_recovered_all = {date_cum_recovered1,date_cum_recovered2};
date_new_confirmed_all = {date_new_confirmed1,date_new_confirmed2};
date_new_recovered_all = {date_new_recovered1,date_new_recovered2};
date_quarantined_all = {date_quarantined1,date_quarantined2};


for i = 1:2
    DataTobeFitted = Data_all{i};
    N_all = NN{i};
    data_pre = data_all{i};
    time_pre = time_all{i};
    date_cum_confirmed = date_cum_confirmed_all{i};
    date_cum_recovered = date_cum_recovered_all{i};
    date_new_confirmed = date_new_confirmed_all{i};
    date_new_recovered = date_new_recovered_all{i};
    date_quarantined   = date_quarantined_all{i};
    
    quarantined_fit{i}   = DataTobeFitted(1:N_all(1));
    cum_confirmed_fit{i} = DataTobeFitted(N_all(1)+1:N_all(1)+N_all(2));
    cum_recovered_fit{i} = DataTobeFitted(N_all(1)+N_all(2)+1:N_all(1)+N_all(2)+N_all(3));
    new_confirmed_fit{i} = DataTobeFitted(N_all(1)+N_all(2)+N_all(3)+1:N_all(1)+N_all(2)+N_all(3)+N_all(4));
    new_recovered_fit{i} = DataTobeFitted(N_all(1)+N_all(2)+N_all(3)+N_all(4)+1:N_all(1)+N_all(2)+N_all(3)+N_all(4)+N_all(5));


    now_qz_pre{i}   = data_pre(:,1);
    cum_qz_pre{i}   = data_pre(:,2);
    cum_cure_pre{i} = data_pre(:,3);
    new_qz_pre{i}   = data_pre(:,4);
    new_cure_pre{i} = data_pre(:,5);
    rate_qz_pre{i}  = data_pre(:,6);
    rate_cure_pre{i} = data_pre(:,7);
    S_pre{i} = data_pre(:,8);
    H_pre{i} = data_pre(:,9);
    
    
%% Caculate the rate

    new_confirmed_fit2 = [0;new_confirmed_fit{i}];
    new_recovered_fit2 = [0;new_recovered_fit{i}];
    
    rate_qz = [];
    rate_cure = [];
    
    for k = 2:N_all(4)+1
        rate_qz(k) = new_confirmed_fit2(k)./quarantined_fit{i}(k-1);
    end

    for k = 2:N_all(5)+1
        rate_cure(k) = new_recovered_fit2(k)./quarantined_fit{i}(k-1);
    end
    
    rate_confirmed{i} = rate_qz;
    rate_recovered{i} = rate_cure;
end


%% disp

disp(['gamma: ' num2str(roundn(gamma(1),-2)),'——',num2str(roundn(gamma(2),-2))])

disp(['alpha: ' num2str(roundn(alpha(1),-2)),'——',num2str(roundn(alpha(2),-2))])

disp(['S1: ' num2str(roundn(s1(1),-2)),'——',num2str(roundn(s1(2),-2))])

disp(['confirmed: ' num2str(cum_cases(1)),'——',num2str(cum_cases(2))])

disp(['The first wave: ' datestr(date_cum_confirmed1(1)) '——' ...
    datestr(date_cum_confirmed1(length(rate_confirmed{1}))) ...
    '(' num2str(length(rate_confirmed{1})) 'days)'])

disp(['The second wave: ' datestr(date_cum_confirmed2(1)) '——' ...
    datestr(date_cum_confirmed2(length(rate_confirmed{2}))) ...
    '(' num2str(length(rate_confirmed{2})) 'days)'])

disp(['The ration of gamma: ' num2str(roundn(gamma(2),-2)/roundn(gamma(1),-2))])
disp(['The ration of alpha: ' num2str(roundn(alpha(2),-2)/roundn(alpha(1),-2))])


%% Plot
gcf = figure('position',[382,107,1186,668]);
t = tiledlayout(10,3,'TileSpacing','tight','Padding','compact');

% Infection process
% Separating degree
gca = nexttile([4,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;

hold on
plot(S_pre{1},'-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(S_pre{2},'-','LineWidth',2,'color',[0.80392	0.33333	0.33333])

xlabel('days','FontName','Times New Roman','FontSize',16)
ylabel('separating degree','FontName','Times New Roman','FontSize',16)

legend({'first wave','second wave'}....
    ,'location','northeast','FontSize',16,'numcolumns',1)
legend('boxoff')

set(gca,'FontSize',16)
box on
xlim([1 50])
ylim([0,1])

title(['(a) '  ' separating degree'],'Position',[25,-0.45,0],'FontWeight','normal',...
    'FontName','Times New Roman','FontSize',18)


% Infection rate
gca = nexttile([4,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;

hold on
plot(rate_confirmed{1}(1:end),'^','MarkerSize',6,'linewidth',1,...
    'color',[0.27451	0.5098	0.70588]);
plot(rate_confirmed{2}(1:end),'^','MarkerSize',6,'linewidth',1,...
    'color',[0.80392	0.33333	0.33333]);

plot(rate_qz_pre{1},'-','LineWidth',2,'color',[0.27451	0.5098	0.70588]);
plot(rate_qz_pre{2},'-','LineWidth',2,'color',[0.80392	0.33333	0.33333]);

xlabel({'days'},...
    'FontName','Times New Roman','FontSize',16)
ylabel('infection rate',...
    'FontName','Times New Roman','FontSize',16)

legend({'first wave','second wave'}...
    ,'location','northeast','FontSize',16,'numcolumns',1)
legend('boxoff')

set(gca,'FontSize',16)
box on
xlim([1 50])
ylim([10^-4 10^2])
yticks([10^-4,10^-2,10^0,10^2])
yticklabels([10^-4,10^-2,10^0,10^2])
ytickformat('%.1e');
set(gca,'yscale','log')

title(['(b) '  ' infection rate'],'Position',[25,0.0000002,0],'FontWeight','normal',...
    'FontName','Times New Roman','FontSize',18)

% confirmed
gca = nexttile([4,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;

hold on;

plot(cum_confirmed_fit{1}(1:length(rate_confirmed{1}(1:end))),'o',...
    'MarkerSize',6,'linewidth',1,'color',[0.27451	0.5098	0.70588]);
plot(cum_confirmed_fit{2}(1:length(rate_confirmed{2}(1:end))),'o',...
    'MarkerSize',6,'linewidth',1,'color',[0.80392	0.33333	0.33333]);

plot(cum_qz_pre{1},'-','MarkerSize',6,'linewidth',2,...
    'color',[0.27451	0.5098	0.70588]);
plot(cum_qz_pre{2},'r-','MarkerSize',6,'linewidth',2,...
    'color',[0.80392	0.33333	0.33333]);

xlabel('days','FontName','Times New Roman','FontSize',16)
ylabel('infected numbers','FontName','Times New Roman','FontSize',16)

legend({'first wave','second wave'},'FontSize',16,'numcolumns',1);
legend('boxoff')
set(gca,'yscale','log')
set(gca,'FontSize',16)
box on
xlim([1 50])
ylim([1,10000])

title(['(c) '  ' infected numbers'],'Position',[25,0.016,0],'FontWeight','normal',...
    'FontName','Times New Roman','FontSize',18)

%%

gca = nexttile([1,3]);
gca.Visible = 'off';

%% Healing process

% Heal force
gca = nexttile([4,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;

hold on

plot(H_pre{1}(gap_new_1:end),'-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(H_pre{2}(gap_new_2:end),'-','LineWidth',2,'color',[0.80392	0.33333	0.33333])

xlabel('days','FontName','Times New Roman','FontSize',16)
ylabel('healing degree','FontName','Times New Roman','FontSize',16)

legend({'first wave','second wave'}....
    ,'location','northeast','FontSize',16,'numcolumns',1)
legend('boxoff')
set(gca,'FontSize',16)
box on
xlim([1 50])
ylim([0,1])

title(['(d) '  ' healing degree'],'Position',[25,-0.45,0],'FontWeight','normal',...
    'FontName','Times New Roman','FontSize',18)


% Cure rate
gca = nexttile([4,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;

hold on
plot(rate_recovered{1}(2:end),'^','MarkerSize',6,'linewidth',1,...
    'color',[0.27451	0.5098	0.70588]);
plot(rate_recovered{2}(2:end),'^','MarkerSize',6,'linewidth',1,...
    'color',[0.80392	0.33333	0.33333]);

plot(rate_cure_pre{1}(gap_new_1:end),'-','LineWidth',2,...
    'color',[0.27451	0.5098	0.70588]);
plot(rate_cure_pre{2}(gap_new_2:end),'-','LineWidth',2,...
    'color',[0.80392	0.33333	0.33333]);

xlabel({'days'},...
    'FontName','Times New Roman','FontSize',16)
ylabel('cure rate',...
    'FontName','Times New Roman','FontSize',16)

legend({'first wave','second wave'}...
    ,'location','northeast','FontSize',16,'numcolumns',1)
legend('boxoff')

set(gca,'FontSize',16)
box on
xlim([1 50])
ylim([10^-4 10^2])
yticks([10^-4,10^-2,10^0,10^2])
yticklabels([10^-4,10^-2,10^0,10^2])
ytickformat('%.1e');
set(gca,'yscale','log')

title(['(e) '  ' cure rate'],'Position',[25,0.0000002,0],'FontWeight','normal',...
    'FontName','Times New Roman','FontSize',18)


% recovered
gca = nexttile([4,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;

hold on;
plot(cum_recovered_fit{1}(1:length(rate_recovered{1}(2:end))),'o',...
    'MarkerSize',6,'linewidth',1,'color',[0.27451	0.5098	0.70588]);
plot(cum_recovered_fit{2}(1:length(rate_recovered{2}(2:end))),'o',...
    'MarkerSize',6,'linewidth',1,'color',[0.80392	0.33333	0.33333]);

plot(cum_cure_pre{1}(gap_cum_1:end),'-','linewidth',2,...
    'color',[0.27451	0.5098	0.70588]);
plot(cum_cure_pre{2}(gap_cum_2:end),'-','linewidth',2,...
    'color',[0.80392	0.33333	0.33333]);

xlabel('days','FontName','Times New Roman','FontSize',16)
ylabel('cured numbers','FontName','Times New Roman','FontSize',16)

L1=legend({'first wave','second wave'},'FontSize',16,'numcolumns',1);
legend('boxoff')
set(gca,'yscale','log')
set(gca,'FontSize',16)
box on
xlim([1 50])
ylim([1,10000])
xlabel('days')


title(['(f) '  ' cured numbers'],'Position',[25,0.016,0],'FontWeight','normal',...
    'FontName','Times New Roman','FontSize',18)

gca = nexttile([1,3]);
gca.Visible = 'off';
