% This file is used to plot the figure 3
% Author: Lei Zhang
% Last modified: 2022-06-27


clear
close all
clc

%% read the data : Hubei
inputData1 = xlsread('hubei.xlsx');
% confirmed
cum_hb_qz = inputData1(:,1);
new_hb_qz = inputData1(:,2);
% deaths
cum_hb_death = inputData1(:,3);
new_hb_death = inputData1(:,4);
% cured
cum_hb_cure = inputData1(:,5);
new_hb_cure = inputData1(:,6);

% quarantined
now_hb = cum_hb_qz - cum_hb_cure - cum_hb_death;

% rate
for i = 2:59 
    rate_hb_qz(i) = new_hb_qz(i)./now_hb(i-1);
end

for i = 2:95 
    rate_hb_cure(i) = new_hb_cure(i)./now_hb(i-1);
end

for i = 2:79
    rate_hb_death(i) = new_hb_death(i)./now_hb(i-1);
end



%% read the data : Italy
inputFilename2 = 'Italy.xlsx';
inputData2 = xlsread(inputFilename2);

% confirmed
cum_italy_qz = inputData2(:,1);
new_italy_qz = inputData2(:,2);
% deaths
cum_italy_death = inputData2(:,3);
new_italy_death = inputData2(:,4);
% cured
cum_italy_cure = inputData2(:,5);
new_italy_cure = inputData2(:,6);

% quarantined
now_italy = cum_italy_qz - cum_italy_cure - cum_italy_death;

% rate
for i = 2:101
    rate_italy_qz(i) = new_italy_qz(i)./now_italy(i-1);
end

for i = 2:101
    rate_italy_cure(i) = new_italy_cure(i)./now_italy(i-1);
end

for i = 2:101
    rate_italy_death(i)  = new_italy_death(i) ./now_italy(i-1);
end



%% the simulation of Hubei
% the initial reported data
I0_hb    = now_hb(1);

% the parameters of separating degree (S)
i0_hb    = 0.48;
gamma_hb = 0.18;
s1_hb    = 1;
ts_hb    = 17.8;

% the parameters of healing degree (H)
alpha_hb = 0.08;
h1_hb    = 0.15;
th_hb    = 50;

% the parameters of rescuing degree (R)
d0_hb   = 0.027;
mu_hb   = 0.16;
r1_hb   = 0.95;
tr_hb   = 10;

% simulation
N_hb = 400;
yhat_hb = OdeNew(N_hb,i0_hb,gamma_hb,s1_hb,alpha_hb,h1_hb,d0_hb,mu_hb,r1_hb,I0_hb,ts_hb,th_hb,tr_hb,...
    cum_hb_qz(1),cum_hb_cure(1),cum_hb_death(1),...
    new_hb_qz(1),new_hb_cure(1),new_hb_death(1));

% the simulation result
now_hbqz_pre =  exp(yhat_hb(:,1));

new_hbqz_pre    =  exp(yhat_hb(:,2));
new_hbcure_pre  =  exp(yhat_hb(:,3));
new_hbdeath_pre =  exp(yhat_hb(:,4));

cum_hbqz_pre    =  exp(yhat_hb(:,5));
cum_hbcure_pre  =  exp(yhat_hb(:,6));
cum_hbdeath_pre =  exp(yhat_hb(:,7));

rate_hbqz_pre    = i0_hb*(1-exp(yhat_hb(:,8)));
rate_hbcure_pre  = exp(yhat_hb(:,9));
rate_hbdeath_pre = d0_hb*(1-exp(yhat_hb(:,10)));

S_hb_pre = exp(yhat_hb(:,8));
H_hb_pre = exp(yhat_hb(:,9));
R_hb_pre = exp(yhat_hb(:,10));

S1_hb  = s1_hb /(1+exp(-gamma_hb *(1-ts_hb )));
H1_hb  = h1_hb/(1+exp(-alpha_hb *(1-th_hb )));
R1_hb  = r1_hb/(1+exp(-mu_hb *(1-tr_hb )));



%% the simulation of Italy
% the initial reported data
I1_italy    = now_italy(1);

% the parameters of separating degree (S)
i0_italy    = 0.48;
gamma_italy = 0.10;
s1_italy    = 0.98;
ts_italy    = 17;

% the parameters of healing degree (H)
alpha_italy = 0.03;
h1_italy    = 0.09;
tz_italy    = 95;

% the parameters of rescuing degree (R)
d0_italy     = 0.02;
mu_italy     = 0.09;
r1_italy     = 0.9;
tr_italy     = 31;

% simulation
N_italy = 400;
yhat_italy = OdeNew(N_italy,i0_italy,gamma_italy,s1_italy,alpha_italy,h1_italy,d0_italy,mu_italy,r1_italy,I1_italy,ts_italy,tz_italy,tr_italy,...
    cum_italy_qz(1),cum_italy_cure(1),cum_italy_death(1),...
    new_italy_qz(1),new_italy_cure(1),new_italy_death(1));

% the simulation result
now_italyqz_pre = exp(yhat_italy(:,1));

new_italyqz_pre = exp(yhat_italy(:,2));
new_italycure_pre = exp(yhat_italy(:,3));
new_italydeath_pre = exp(yhat_italy(:,4));

cum_italyqz_pre    =  exp(yhat_italy(:,5));
cum_italycure_pre  =  exp(yhat_italy(:,6));
cum_italydeath_pre =  exp(yhat_italy(:,7));

rate_italyqz_pre    = i0_italy*(1-exp(yhat_italy(:,8)));
rate_italycure_pre  =  exp(yhat_italy(:,9));
rate_italydeath_pre =  d0_italy*(1-exp(yhat_italy(:,10)));

S_italy_pre = exp(yhat_italy(:,8));
H_italy_pre = exp(yhat_italy(:,9));
R_italy_pre = exp(yhat_italy(:,10));

S1_italy  = s1_italy/(1+exp(-gamma_italy *(1-ts_italy )));
H1_italy  = h1_italy/(1+exp(-alpha_italy *(1-tz_italy )));
R1_italy  = r1_italy/(1+exp(-mu_italy *(1-tr_italy )));

%% Plot the result

figure('position',[382,107,1186,809]);
t = tiledlayout(22,3,'TileSpacing','tight','Padding','compact');

% the parameters of separating degree (S)
gca = nexttile([6,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on
plot(S_hb_pre,'-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(S_italy_pre,'-','LineWidth',2,'color',[0.80392	0.33333	0.33333])

xlabel('days','FontName','Times New Roman','FontSize',16)
ylabel('separating degree','FontName','Times New Roman','FontSize',16)
legend({'Hubei','Italy'},'FontSize',16,...
    'location','east','NumColumns',1,'interpreter','latex')
legend('boxoff')
set(gca,'FontSize',16)
box on
set(gca,'xtick',[1  50 100 150]);
set(gca,'xticklabel',{'1','50','100','150'});
axis([1 150 0 1])
title(['$(a)$ $\quad$ separating degree $(S)$'],'Position',[75,-0.65,0],...
    'FontSize',18,'interpreter','latex')


% the parameters of healing degree (H)
gca = nexttile([6,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on
plot(H_hb_pre,'-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(H_italy_pre,'-','LineWidth',2,'color',[0.80392	0.33333	0.33333])

xlabel({'days'},'FontName','Times New Roman','FontSize',16)
ylabel('healing degree','FontName','Times New Roman','FontSize',16)

legend({'Hubei','Italy'},'FontSize',16,...
    'location','east','NumColumns',1,'interpreter','latex')
legend('boxoff')
set(gca,'FontSize',16)
box on
set(gca,'xtick',[1  50 100 150]);
set(gca,'xticklabel',{'1','50','100','150'});
axis([1 150 0 1])
title(['$(b)$ $\quad$ healing degree $(H)$'],'Position',[75,-0.65,0],...
    'FontSize',18,'interpreter','latex')

% the parameters of rescuing degree (R)
gca = nexttile([6,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on
plot(R_hb_pre,'-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(R_italy_pre,'-','LineWidth',2,'color',[0.80392	0.33333	0.33333])

xlabel('days','FontName','Times New Roman','FontSize',16)
ylabel('rescuing degree','FontName','Times New Roman','FontSize',16)

legend({'Hubei','Italy'},'FontSize',16,...
    'location','east','NumColumns',1,'interpreter','latex')
legend('boxoff')
set(gca,'FontSize',16)
box on
set(gca,'xtick',[1  50 100 150]);
set(gca,'xticklabel',{'1','50','100','150'});
axis([1 150 0 1])

title(['$(c)$ $\quad$ rescuing degree $(R)$'],'Position',[75,-0.65,0],...
    'FontSize',18,'interpreter','latex')

%%

gca = nexttile([1,3]);
gca.Visible = 'off';

%% the transition rate

% infection rate
gca = nexttile([6,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on;
plot(rate_hb_qz,'^','MarkerSize',6,'color',[0.27451	0.5098	0.70588]);
plot(rate_italy_qz,'^','MarkerSize',6,'color',[0.80392	0.33333	0.33333]);
plot(rate_hbqz_pre,'-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(rate_italyqz_pre,'-','LineWidth',2,'color',[0.80392	0.33333	0.33333])

xlabel('days','FontName','Times New Roman','FontSize',16)
ylabel('infection rate','FontSize',16)

legend({'Hubei','Italy'},'FontSize',16,...
    'location','northeast','NumColumns',1,'interpreter','latex')
legend('boxoff')
set(gca,'yscale','log')
set(gca,'FontSize',16)
box on
set(gca,'xtick',[1  50 100 150]);
set(gca,'xticklabel',{'1','50','100','150'});
set(gca,'ytick',[0.00001 0.001  0.1  10 100 1000 10000 100000 1000000]);
axis([1 150 0.00001 50])

title(['$(d)$ $\quad$ infection rate $(V_I)$'],'Position',[75,5e-10,0],...
    'FontSize',18,'interpreter','latex')

% cure rate
gca = nexttile([6,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on
plot(rate_hb_cure,'^','MarkerSize',6,'color',[0.27451	0.5098	0.70588]);
plot(rate_italy_cure,'^','MarkerSize',6,'color',[0.80392	0.33333	0.33333]);
plot(rate_hbcure_pre,'-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(rate_italycure_pre,'-','LineWidth',2,'color',[0.80392	0.33333	0.33333])

xlabel({'days'},'FontName','Times New Roman','FontSize',16)
ylabel('cure rate','FontName','Times New Roman','FontSize',16)

legend({'Hubei','Italy'},'FontSize',16,...
    'location','southeast','NumColumns',1,'interpreter','latex')
legend('boxoff')
set(gca,'yscale','log')
set(gca,'FontSize',16)
box on
set(gca,'xtick',[1  50 100 150]);
set(gca,'xticklabel',{'1','50','100','150'});
set(gca,'ytick',[0.0001  0.01  1 10 100 1000 10000 100000 1000000]);
axis([1 150 0.0001 1])


title(['$(e)$ $\quad$ cure rate $(V_C)$'],'Position',[75,3e-7,0],...
    'FontSize',18,'interpreter','latex')

% death rate
gca = nexttile([6,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on
plot(rate_hb_death,'^','MarkerSize',6,'color',[0.27451	0.5098	0.70588]);
plot(rate_italy_death,'^','MarkerSize',6,'color',[0.80392	0.33333	0.33333]);
plot(rate_hbdeath_pre,'-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(rate_italydeath_pre,'-','LineWidth',2,'color',[0.80392	0.33333	0.33333])

xlabel('days','FontName','Times New Roman','FontSize',16)
ylabel('death rate','FontName','Times New Roman','FontSize',16)

legend({'Hubei','Italy'},'FontSize',16,...
    'location','northeast','interpreter','latex')
legend('boxoff')
set(gca,'yscale','log')
set(gca,'FontSize',16)
box on
set(gca,'xtick',[1  50 100 150]);
set(gca,'xticklabel',{'1','50','100','150'});
set(gca,'ytick',[0.0001  0.001 0.01 0.1  1 10 100 1000 10000 100000 1000000]);
axis([1 150 0.0001 0.3])


title(['$(f)$ $\quad$ death rate $(V_D)$'],'Position',[75,6e-7,0],...
    'FontSize',18,'interpreter','latex')



%% the report cases

gca = nexttile([1,3]);
gca.Visible = 'off';

%% new cases
time_pre = datetime('2020-01-18'):1:datetime('2020-01-18')+500;
% new confirmed
gca = nexttile([6,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on
plot(time_pre(1:length(new_hb_qz)),new_hb_qz,...
    'o','MarkerSize',6,'color',[0.27451	0.5098	0.70588]);
plot(time_pre(36:1:length(new_italy_qz)+35),new_italy_qz,...
    'o','MarkerSize',6,'color',[0.80392	0.33333	0.33333]);
plot(time_pre(1:length(new_hbqz_pre)),new_hbqz_pre,...
    '-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(time_pre(36:1:length(new_italyqz_pre)+35),new_italyqz_pre,...
    '-','LineWidth',2,'color',[0.80392	0.33333	0.33333])

% xlabel('date','FontName','Times New Roman','FontSize',16)
ylabel('daily infected','FontName','Times New Roman','FontSize',16)

legend({'Hubei','Italy'},'FontSize',16,...
    'location','northeast','NumColumns',1 ,'interpreter','latex')
legend('boxoff')
set(gca,'yscale','log')
set(gca,'FontSize',16)
box on

set(gca,'ytick',[0.0001 0.001 0.01 0.1 1  100  10000 1000000 100000000]);
xlim([time_pre(2),time_pre(152)])
ylim([1,1000000])

title(['$(g)$ $\quad$ daily infected $(\Delta I)$'],'Position',[75,5e-5,0],...
    'FontSize',18,'interpreter','latex')



% new cured
gca = nexttile([6,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on
plot(time_pre(1:length(new_hb_cure)),new_hb_cure,...
    'o','MarkerSize',6,'color',[0.27451	0.5098	0.70588]);
plot(time_pre(36:1:length(new_italy_cure)+35),new_italy_cure,...
    'o','MarkerSize',6,'color',[0.80392	0.33333	0.33333]);
plot(time_pre(1:length(new_hbcure_pre)),new_hbcure_pre,...
    '-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(time_pre(36:1:length(new_italycure_pre)+35),new_italycure_pre,...
    '-','LineWidth',2,'color',[0.80392	0.33333	0.33333])

% xlabel({'date'},'FontName','Times New Roman','FontSize',16)
ylabel('daily cured','FontName','Times New Roman','FontSize',16)

legend({'Hubei','Italy'},'FontSize',16,...
    'location','northeast','NumColumns',1 ,'interpreter','latex')
legend('boxoff')
set(gca,'yscale','log')
set(gca,'FontSize',16)
box on

set(gca,'ytick',[0.0001 0.001 0.01 0.1 1  100  10000 1000000 100000000]);
xlim([time_pre(2),time_pre(152)])
ylim([1,1000000])


title(['$(h)$ $\quad$ daily cured $(\Delta C)$'],'Position',[75,5e-5,0],...
    'FontSize',18,'interpreter','latex')



% new deaths
gca = nexttile([6,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on
plot(time_pre(1:length(new_hb_death)),new_hb_death,...
    'o','MarkerSize',6,'color',[0.27451	0.5098	0.70588]);
plot(time_pre(36:1:length(new_italy_death)+35),new_italy_death,...
    'o','MarkerSize',6,'color',[0.80392	0.33333	0.33333]);
plot(time_pre(1:length(new_hbdeath_pre)),new_hbdeath_pre,...
    '-','LineWidth',2,'color',[0.27451	0.5098	0.70588])
plot(time_pre(36:1:length(new_italydeath_pre)+35),new_italydeath_pre,'-',...
    'LineWidth',2,'color',[0.80392	0.33333	0.33333])

% xlabel('date','FontName','Times New Roman','FontSize',16)
ylabel('daily deaths','FontName','Times New Roman','FontSize',16)

legend({'Hubei','Italy'},'FontSize',16,...
    'location','northeast','NumColumns',1 ,'interpreter','latex')
legend('boxoff')
set(gca,'yscale','log')
set(gca,'FontSize',16)
box on

set(gca,'ytick',[0.0001 0.001 0.01 0.1 1  100  10000 1000000 100000000]);
xlim([time_pre(2),time_pre(152)])
ylim([1,1000000])
title(['$(i)$ $\quad$ daily deaths $(\Delta D)$'],'Position',[75,5e-5,0],...
    'FontSize',18,'interpreter','latex')




%% the model
function yhat = OdeNew(N, i0, gamma,s1,alpha, h1,d0,mu,r1,I1,ts,th,tr,...
    cum_I1,cum_C1,cum_D1,new_I1,new_C1,new_D1)


I_pre=zeros(N,1);
S_pre=zeros(N,1);
H_pre=zeros(N,1);
R_pre=zeros(N,1);
dI_pre=zeros(N,1);
dC_pre=zeros(N,1);
dD_pre=zeros(N,1);
cum_I_pre=zeros(N,1);
cum_C_pre=zeros(N,1);
cum_D_pre=zeros(N,1);

I_pre(1)=I1;

S_1 = s1/(1+exp(-gamma*(1-ts)));
H_1 = h1/(1+exp(-alpha*(1-th)));
R_1 = r1/(1+exp(-mu*(1-tr)));

S_pre(1)=S_1;
H_pre(1)=H_1;
R_pre(1)=R_1;

dI_pre(1) = new_I1;
dC_pre(1) = new_C1;
dD_pre(1) = new_D1;

cum_I_pre(1)=cum_I1;
cum_C_pre(1)=cum_C1;
cum_D_pre(1)=cum_D1;

for i=2:1:N
    S_pre(i) = s1/(1+exp(-gamma*(i-ts)));
    H_pre(i) = h1/(1+exp(-alpha*(i-th)));
    R_pre(i) = r1/(1+exp(-mu*(i-tr)));
    
    dI_pre(i)=I_pre(i-1)*i0*(1-S_pre(i));
    dC_pre(i)=I_pre(i-1)*H_pre(i);
    dD_pre(i)=I_pre(i-1)*d0*(1-R_pre(i));
    
    cum_I_pre(i)=cum_I_pre(i-1)+dI_pre(i);
    cum_C_pre(i)=cum_C_pre(i-1)+dC_pre(i);
    cum_D_pre(i)=cum_D_pre(i-1)+dD_pre(i);
    
    I_pre(i)=I_pre(i-1)+dI_pre(i)-dC_pre(i)-dD_pre(i);
end
    yhat=log([I_pre,dI_pre,dC_pre,dD_pre,cum_I_pre,cum_C_pre,cum_D_pre,S_pre,H_pre,R_pre]);
end