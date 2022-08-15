% This file is used to plot the figure 7
% Author: Lei Zhang
% Last modified: 2022-06-27

clear
close all
clc
inputarea = 'Beijing';

% read the model parameters
tableBJ = readtable('Beijing_diffDataset.csv');
s1 = tableBJ.s1;
cum_I_pre = tableBJ.confirmed_pre;


% the parameters bound of separating degree S(i0, gamma, ts, s1)
% loose interventin: s1 = 93%
para_lb1 = [0,  0,   0,      0.9300];
para_ub1 = [3,  1,   0.001,  0.9300];

% current: free bound
para_lb2 = [0,  0,  0,       0];
para_ub2 = [3,  1,  0.001,   1];

% tight interventin: s1 = 98%
para_lb3 = [0,  0,  0,      0.9795];
para_ub3 = [3,  1,  0.001,  0.9795];

% Call the simulation function
% simulate the wave
num_wave = 3;

[Data1,data_pre1,N_1,time_pre1,date_cum_confirmed1,date_cum_recovered1,...
    date_new_confirmed1,date_new_recovered1,date_quarantined1,...
    para_out1] = func_simulate_wave(num_wave,inputarea,para_lb1,para_ub1);

[Data2,data_pre2,N_2,time_pre2,date_cum_confirmed2,date_cum_recovered2,...
    date_new_confirmed2,date_new_recovered2,date_quarantined2,...
    para_out2] = func_simulate_wave(num_wave,inputarea,para_lb2,para_ub2);

[Data3,data_pre3,N_3,time_pre3,date_cum_confirmed3,date_cum_recovered3,...
    date_new_confirmed3,date_new_recovered3,date_quarantined3,...
    para_out3] = func_simulate_wave(num_wave,inputarea,para_lb3,para_ub3);

% fit with different dataset sizes
% as of May 12
[Data4,data_pre4,N_4,time_pre4,date_cum_confirmed4,date_cum_recovered4,...
    date_new_confirmed4,date_new_recovered4,date_quarantined4]...
    = func_fit_3wave(12,inputarea);

% as of May 21
[Data5,data_pre5,N_5,time_pre5,date_cum_confirmed5,date_cum_recovered5,...
    date_new_confirmed5,date_new_recovered5,date_quarantined5]...
    = func_fit_3wave(3,inputarea);


Data_all = {Data1,Data2,Data3,Data4,Data5};
NN = {N_1,N_2,N_3,N_4,N_5};
data_all = {data_pre1,data_pre2,data_pre3,data_pre4,data_pre5};
time_all = {time_pre1,time_pre2,time_pre3,time_pre4,time_pre5};

date_cum_confirmed_all = {date_cum_confirmed1,date_cum_confirmed2,...
    date_cum_confirmed3,date_cum_confirmed4,date_cum_confirmed5};
date_cum_recovered_all = {date_cum_recovered1,date_cum_recovered2,...
    date_cum_recovered3,date_cum_recovered4,date_cum_recovered5};

date_new_confirmed_all = {date_new_confirmed1,date_new_confirmed2,...
    date_new_confirmed3,date_new_confirmed4,date_new_confirmed5};
date_new_recovered_all = {date_new_recovered1,date_new_recovered2,...
    date_new_recovered3,date_new_recovered4,date_new_recovered5};

date_quarantined_all = {date_quarantined1,date_quarantined2,...
    date_quarantined3,date_quarantined4,date_quarantined5};


num_case(1) = 0;
    
for i = 1:5
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
    
    
    % predict the end date and total confirmed
    date_pre = date_cum_confirmed(1):1:date_cum_confirmed(1)+100;
    [~,ind_newmax] = max(new_qz_pre{i});
    ind_newend1 = find(new_qz_pre{1,i}(ind_newmax:end)<10);
    ind_newend2 = ind_newend1(1) + ind_newmax -1;
    date_newend1 = date_pre(ind_newend2); 

    ind_newend3 = find(new_qz_pre{1,i}(ind_newmax:end)<1);
    ind_newend4 = ind_newend3(1) + ind_newmax -1;
    
    num_case(i+1)=round(cum_qz_pre{1,i}(ind_newend2));
    
    disp(['less than 10: ' datestr(date_newend1)])
    disp(['total confirmed：' num2str(num_case(i+1))])
    disp(['add：' num2str(num_case(i+1)-num_case(i))])
    
end

disp(['average：' num2str(mean(cum_I_pre)),'±', num2str(std(cum_I_pre))])



%% calculate the MRE
data_pre2 = [1404;1441;1464;1484;1498;1506;1530;1550;1564;1572;1588;1602;...
    1613;1621;1630;1647;1652;1656];

realData = [cum_confirmed_fit{1,3};data_pre2];
predResult =  cum_qz_pre{1,5}(1:length(realData));
mre = abs(predResult-realData)./realData;
mre_mean = mean(mre);

disp(['MRE(as of 5.21)：' num2str(roundn(mre_mean,-4)*100) '%'])

%% Plot
gcf = figure('position',[382,107,1186,809]);
t = tiledlayout(12,6,'TileSpacing','compact','Padding','compact');

% xtick
xticks1 = time_pre(1):45:time_pre(1)+90;

% separating degree
gca = nexttile([5,2]);
gca.FontName ='Times New Roman';
gca.FontSize =12;

hold on
plot(time_pre,S_pre{1},'-','LineWidth',2,'color',[0.80392	0.33333	0.33333])
plot(time_pre,S_pre{2},'-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(time_pre,S_pre{3},'-','LineWidth',2,'color',[0.47451	0.80392	0.80392])

ylabel('separating degree','FontName','Times New Roman','FontSize',18)

legend({['loose      S_1 = ',num2str(roundn(para_out1(5),-3)*100) '%'],....
    ['23 May  S_1 = ',num2str(roundn(para_out2(5),-3)*100) '%'],...
    ['tight       S_1 = ',num2str(roundn(para_out3(5),-3)*100) '%']},...
    'location','southeast','FontSize',14,'numcolumns',1)
legend('boxoff')
box on
xlim([time_pre(1) xticks1(end)])
ylim([0.5,1])
set(gca,'FontSize',16)

title(['(a) '  ' separating degree'],...
    'Position',[45,0.33,0],'FontWeight','normal',...
    'FontName','Times New Roman','FontSize',18);

%% Infection rate

gca = nexttile([5,2]);
gca.FontName ='Times New Roman';
gca.FontSize =12;

hold on
plot(date_new_confirmed_all{1},rate_confirmed{1}(2:end),'k^','MarkerSize',8,'linewidth',1);
plot(time_pre,rate_qz_pre{1},'-','LineWidth',2,'color',[0.80392	0.33333	0.33333]);
plot(time_pre,rate_qz_pre{2},'-','LineWidth',2,'color',[0.27451	0.5098	0.70588]);
plot(time_pre,rate_qz_pre{3},'-','LineWidth',2,'color',[0.47451	0.80392	0.80392]);


ylabel('infection rate',...
    'FontName','Times New Roman','FontSize',18)

legend({'data as of 23 May' ['loose      S_1 = ',num2str(roundn(para_out1(5),-3)*100) '%'],....
    ['23 May  S_1 = ',num2str(roundn(para_out2(5),-3)*100) '%'],...
    ['tight       S_1 = ',num2str(roundn(para_out3(5),-3)*100) '%']},...
    'location','northeast','FontSize',14,'numcolumns',1);
legend('boxoff')

set(gca,'yscale','log')
box on
xlim([time_pre(1) xticks1(end)])
ylim([10^-2 10^1])
yticks([10^-2,10^-1,10^0,10^1])
yticklabels([10^-2,10^-1,10^0,10^1])
ytickformat('%.1e');
set(gca,'FontSize',16)


title(['(b) '  ' infection rate'],...
    'Position',[45,0.0009,0],'FontWeight','normal',...
    'FontName','Times New Roman','FontSize',18);

%% Total confirmed
gca = nexttile([5,2]);
gca.FontName ='Times New Roman';
gca.FontSize =12;

hold on;
plot(date_cum_confirmed_all{1},cum_confirmed_fit{1},'ko','MarkerSize',8,'linewidth',1);
plot(time_pre,cum_qz_pre{1},'c-','linewidth',2,'color',[0.80392	0.33333	0.33333]);
plot(time_pre,cum_qz_pre{2},'b-','linewidth',2,'color',[0.27451	0.5098	0.70588]);
plot(time_pre,cum_qz_pre{3},'r-','linewidth',2,'color',[0.47451	0.80392	0.80392]);

ylabel('total infected','FontName','Times New Roman','FontSize',18)

legend({'data as of 23 May' ['loose      S_1 = ',num2str(roundn(para_out1(5),-3)*100) '%'],....
    ['23 May  S_1 = ',num2str(roundn(para_out2(5),-3)*100) '%'],...
    ['tight       S_1 = ',num2str(roundn(para_out3(5),-3)*100) '%']},...
    'location','northwest','FontSize',14,'numcolumns',1);
legend('boxoff')

box on
xlim([time_pre(1) xticks1(end)])
ylim([1,5000])
set(gca,'FontSize',16)

title(['(c) '  ' total infected'],...
    'Position',[45,-1700,0],'FontWeight','normal',...
    'FontName','Times New Roman','FontSize',18);

%%

gca = nexttile([1,6]);
gca.Visible = 'off';


%%
gca = nexttile([5,3]);
gca.FontName ='Times New Roman';
gca.FontSize =12;

hold on;
plot(cum_I_pre,s1,'o-','MarkerSize',8,'linewidth',2,...
    'color',[0.27451	0.5098	0.70588],'Markerfacecolor',[0.27451	0.5098	0.70588]);

xlabel('predicted total infected','FontName','Times New Roman','FontSize',18)
ylabel('separating degree','FontName','Times New Roman','FontSize',18)

set(gca,'FontSize',16)
box on

xlim([1300,2200])
ylim([0.963,0.9675])

box on
set(gca,'ytick',[0.963,0.965,0.967,0.9675])

yticklabels({'96.3%','96.5%','96.7%'})


annotation('arrow',[0.3820,0.3050],[0.2800,0.3320],'color',[0.27451	0.5098	0.70588]);
annotation('arrow',[0.2600,0.2070],[0.3655,0.3810],'color',[0.27451	0.5098	0.70588]);
annotation('arrow',[0.2500,0.2685],[0.2450,0.2320],'color',[0.27451	0.5098	0.70588]);

text_dateG = {...
    'May 12','May 13','May 14','May 15','May 16','May 17',...
    'May 18','May 19','May 20','May 21','May 22','May 23'};

text_xy = [...
0.965331860408163	2000
0.966365097591837	1727
0.966657887653061	1561
0.966777662459535   1369
0.966132338216748	1340
0.965565885122449	1347
0.965280862163265	1388
0.964991451979592	1436
0.964690898102041	1485
0.964348896326531	1569
0.963793794285715	1718
0.963774090510204	1913];

for te = 1:12
    text(text_xy(te,2),text_xy(te,1),text_dateG{te},...
        'LineWidth',0.5,'FontSize',10,'color',[0.6 0.6 0.6])
end


title(['(d) '  ' phase trajectory of saturated isolation'],...
    'Position',[1750,0.9614,0],'FontWeight','normal',...
    'FontName','Times New Roman','FontSize',18);

%% 
gca = nexttile([5,3]);
gca.FontName ='Times New Roman';
gca.FontSize =12;

hold on;

date_pre2 = date_cum_confirmed1(end)+1:1:date_cum_confirmed1(end)+length(data_pre2);


plot(date_cum_confirmed_all{1},cum_confirmed_fit{1},'ko','MarkerSize',8,...
    'linewidth',1);

plot(time_pre,cum_qz_pre{4},'-','linewidth',2,'color',[0.95686	0.64314	0.37647])
plot(time_pre,cum_qz_pre{2},'-','linewidth',2,'color',[0.27451	0.5098	0.70588]);
plot(time_pre,cum_qz_pre{5},'-','linewidth',2,'color',[0.57647	0.43922	0.85882]);
plot(date_pre2,data_pre2,'ko','MarkerSize',8,'linewidth',1);

% xlabel('date','FontName','Times New Roman','FontSize',18)
ylabel('total infected','FontSize',18)

legend({'data as of 10 June',...
    'S_1 = 96.5%(as of May 12)',...
    'S_1 = 96.4%(as of May 23)',...
    'S_1 = 96.5%(as of May 21)'},...
    'location','southeast','FontSize',14,'numcolumns',1);

legend('boxoff')

box on
xlim([time_pre(1) xticks1(end)])
ylim([1,2500])
set(gca,'FontSize',16)

title(['(e) '  ' predicted total infected'],...
    'Position',[45,-850,0],'FontWeight','normal',...
    'FontName','Times New Roman','FontSize',18)
