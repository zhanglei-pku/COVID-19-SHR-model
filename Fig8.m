
clear
close all
clc



%% read the data file
inputFilename_zhejiang = 'zhejiang.xlsx';
inputFileName_zhejiang_Compar = 'zhejiang_Compare.xlsx';
inputFilename_hubei = 'hubei.xlsx';
inputFileName_buhei_Compar1 = 'hubei_Compare1.xlsx';
inputFileName_buhei_Compar2 = 'hubei_Compare2.xlsx';

inputData_zhejiang = xlsread(inputFilename_zhejiang);
inputData_zhejiang_Compar = xlsread(inputFileName_zhejiang_Compar);

inputData_hubei = xlsread(inputFilename_hubei);
inputData_buhei_Compar1 = xlsread(inputFileName_buhei_Compar1);
inputData_buhei_Compar2 = xlsread(inputFileName_buhei_Compar2);


now_infect_zhejiang_Compar = smooth(inputData_zhejiang_Compar(:,2),3);
now_infect_buhei_Compar1 = smooth(inputData_buhei_Compar1(:,2),3);
now_infect_buhei_Compar2 = smooth(inputData_buhei_Compar2(:,2),3);

%% Zhejiang province
% infected cases
cum_infect_zhejiang = inputData_zhejiang(:,1);
new_infect_zhejiang = inputData_zhejiang(:,2);

% deaths cases
cum_death_zhejiang = inputData_zhejiang(:,3);
new_death_zhejiang = inputData_zhejiang(:,4);

% cured cases
cum_cure_zhejiang = inputData_zhejiang(:,5);
new_cure_zhejiang = inputData_zhejiang(:,6);

% the imported cases from abroad
cum_infect_zhejiang_import = inputData_zhejiang(:,7);
cum_cure_zhejiang_import = inputData_zhejiang(:,8);
cum_death_zhejiang_import = inputData_zhejiang(:,9);

% quarantined cases still in hospitals
% note that remove the imported cases from abroad
now_infect_zhejiang = cum_infect_zhejiang - cum_cure_zhejiang - cum_death_zhejiang...
    -(cum_infect_zhejiang_import - cum_cure_zhejiang_import - cum_death_zhejiang_import);


% simulation by model
% separating degree
Vi0_zhejiang = 0.56;
gamma_zhejiang = 0.30;
ts_zhejiang = 5.8;

% healing degree
alpha_zhejiang = 0.11;
h1_zhejiang = 0.43;
th_zhejiang = 41;

% time
N_zhejiang = 400;

% initial value
% quarantined
Q1_zhejiang    = 102;

cum_I1_zhejiang = cum_infect_zhejiang(1);
cum_C1_zhejiang = cum_cure_zhejiang(1);

new_I1_zhejiang = new_infect_zhejiang(1);
new_C1_zhejiang = new_cure_zhejiang(1);

yhat = OdeNew_zhejiang(N_zhejiang,Vi0_zhejiang,gamma_zhejiang,...
    alpha_zhejiang,h1_zhejiang,Q1_zhejiang,ts_zhejiang,th_zhejiang,...
    cum_I1_zhejiang,cum_C1_zhejiang,new_I1_zhejiang,new_C1_zhejiang);

now_infect_zhejiang_pre = exp(yhat(1:400,:));

time_pre_zhejiang = datetime('2020-01-26'):1:datetime('2020-01-26')+200;
time_pre_zhejiang_Compar1= datetime('2020-01-26')+inputData_zhejiang_Compar(:,1);

%% Hubei province
cum_infect_hubei = inputData_hubei(:,1);
new_infect_hubei = inputData_hubei(:,2);

cum_death_hubei = inputData_hubei(:,3);
new_death_hubei = inputData_hubei(:,4);

cum_cure_hubei = inputData_hubei(:,5);
new_cure_hubei = inputData_hubei(:,6);

now_infect_noImport_hubei = cum_infect_hubei-cum_cure_hubei-cum_death_hubei;

Q1_hubei    = 94;
Vi0_hubei = 0.48;
gamma_hubei = 0.18;
s1_hubei = 1;
ts_hubei = 17.8;


alpha_hubei = 0.08;
h1_hubei = 0.15;
th_hubei = 50;


VD0_hubei = 0.027;
mu_hubei = 0.16;
R1_hubei = 0.95;
tr_hubei = 10;


N_hubei = 300;


cum_I0_hubei = cum_infect_hubei(1);
cum_C0_hubei = cum_cure_hubei(1);
cum_D0_hubei = cum_death_hubei(1);
new_I0_hubei = new_infect_hubei(1);
new_C0_hubei = new_cure_hubei(1);
new_D0_hubei = new_death_hubei(1);

yhat = OdeNew_hubei(N_hubei,Vi0_hubei,gamma_hubei,s1_hubei,...
    alpha_hubei,h1_hubei,VD0_hubei,mu_hubei,R1_hubei,Q1_hubei,...
    ts_hubei,th_hubei,tr_hubei,cum_I0_hubei,cum_C0_hubei,cum_D0_hubei,...
    new_I0_hubei,new_C0_hubei,new_D0_hubei);


now_infect_hubei_pre    =  exp(yhat(1:400,:));

time_pre_hubei=datetime('2020-01-19'):1:datetime('2020-01-19')+200;

time_pre_hubei_Compar1= datetime('2020-01-19')+inputData_buhei_Compar1(:,1);
time_pre_hubei_Compar2= datetime('2020-01-19')+inputData_buhei_Compar2(:,1);

%% plot
gcf = figure('position',[100,100,1186,469]);
t = tiledlayout(10,2,'TileSpacing','compact','Padding','compact');

gca = nexttile([8,1]);
gca.FontName = 'Times New Roman';
gca.FontSize = 12;
hold on


plot(time_pre_zhejiang(1:length(now_infect_zhejiang)),...
    now_infect_zhejiang,...
    'ko','MarkerSize',8,'LineWidth',0.5);
plot(time_pre_zhejiang(1:100),...
    now_infect_zhejiang_pre(1:100),...
    '-','LineWidth',3,'color',[0.80392	0.33333	0.33333])
plot(time_pre_zhejiang_Compar1(1:end)-5,...
    now_infect_zhejiang_Compar(1:end),...
    '-','LineWidth',3,'color',[0.27451	0.5098	0.70588]);

xlim([time_pre_zhejiang(1),time_pre_zhejiang(53)])
ylim([1,1e4])

ylabel('quarantined  cases')
legend({'data','SHR','SEIR$^{[13]}$'},...
    'location','northeast','NumColumns',1,'interpreter','latex')
legend('boxoff')

set(gca,'yscale','log')
set(gca,'FontSize',16)
box on
title(['$(a)$ \quad Zhejiang province'],'Position',[25,0.065,0],...
    'FontSize',18,'interpreter','latex')






gca = nexttile([8,1]);
gca.FontName ='Times New Roman';
gca.FontSize =12;
hold on

plot(time_pre_hubei(1:length(now_infect_noImport_hubei)),...
    now_infect_noImport_hubei,'ko','MarkerSize',8,'LineWidth',0.5);

plot(time_pre_hubei(1:100),now_infect_hubei_pre(1:100),...
    '-','LineWidth',3,'color',[0.80392	0.33333	0.33333])
plot(time_pre_hubei_Compar1(1:42)-6,...
    now_infect_buhei_Compar1(1:42),...
'-','LineWidth',3,'color',[0.27451	0.5098	0.70588]);
plot(time_pre_hubei_Compar2(1:end)+2,...
    now_infect_buhei_Compar2(1:end),...
    '-','LineWidth',3,'color',[0.95686	0.64314	0.37647]);

xlim([time_pre_hubei(1),time_pre_hubei(100)])
ylim([100,1e5])

ylabel('quarantined  cases')
legend({'data','SHR','SEIR$^{[13]}$','SEIR$^{[14]}$'},...
    'location','northeast','interpreter','latex')
legend('boxoff')
set(gca,'yscale','log')
set(gca,'FontSize',16)
box on

title(['$(b)$ \quad Hubei province'],'Position',[45,13,0],...
    'FontSize',18,'interpreter','latex')


%% model
function yhat=OdeNew_zhejiang(N,c0,gamma,alpha,z0,I1,tc,tz,cum_I1,cum_R1,new_I1,new_R1)

I_pre=zeros(N,1);
G_pre=zeros(N,1);
Z_pre=zeros(N,1);
dI_pre=zeros(N,1);
dR_pre=zeros(N,1);
cum_I_pre=zeros(N,1);
cum_R_pre=zeros(N,1);

G1 = 1/(1+exp(-gamma*(1-tc)));
Z1 = z0/(1+exp(-alpha*(1-tz)));

I_pre(1)=I1;
G_pre(1)=G1;
Z_pre(1)=Z1;
dI_pre(1)    = new_I1;
dR_pre(1)    = new_R1;
cum_I_pre(1) = cum_I1;
cum_R_pre(1) = cum_R1;

for i=2:1:400
    G_pre(i) = 1/(1+exp(-gamma*(i-tc)));
    Z_pre(i) = z0/(1+exp(-alpha*(i-tz)));
    
    dI_pre(i)=I_pre(i-1)*c0*(1-G_pre(i));
    dR_pre(i)=I_pre(i-1)*Z_pre(i);
    
    cum_I_pre(i)=cum_I_pre(i-1)+dI_pre(i);
    cum_R_pre(i)=cum_R_pre(i-1)+dR_pre(i);
    
    I_pre(i)=I_pre(i-1)+dI_pre(i)-dR_pre(i);
end

C_pre = c0.*(1-G_pre);

yhat=log([I_pre;cum_I_pre;cum_R_pre;...
    dI_pre;dR_pre;C_pre;Z_pre]);
end

function yhat=OdeNew_hubei(N,c0,gamma,g0,alpha,z0,s0,mu,k0,I1,tc,tz,ts,cum_I1,cum_R1,cum_D1,new_I1,new_R1,new_D1)

I_pre=zeros(N,1);
G_pre=zeros(N,1);
Z_pre=zeros(N,1);
K_pre=zeros(N,1);
dI_pre=zeros(N,1);
dR_pre=zeros(N,1);
dD_pre=zeros(N,1);
cum_I_pre=zeros(N,1);
cum_R_pre=zeros(N,1);
cum_D_pre=zeros(N,1);

I_pre(1)=I1;
G1 = g0/(1+exp(-gamma*(1-tc)));
Z1 = z0/(1+exp(-alpha*(1-tz)));
K1 = k0/(1+exp(-mu*(1-ts)));

G_pre(1)=G1;
Z_pre(1)=Z1;
K_pre(1)=K1;
dI_pre(1)    = new_I1;
dR_pre(1)    = new_R1;
dD_pre(1)    = new_D1;
cum_I_pre(1) = cum_I1;
cum_R_pre(1) = cum_R1;
cum_D_pre(1) = cum_D1;

for i=2:1:400
    G_pre(i) = g0/(1+exp(-gamma*(i-tc)));
    Z_pre(i) = z0/(1+exp(-alpha*(i-tz)));
    K_pre(i) = k0/(1+exp(-mu*(i-ts)));
    
    dI_pre(i)=I_pre(i-1)*c0*(1-G_pre(i));
    dR_pre(i)=I_pre(i-1)*Z_pre(i);
    dD_pre(i)=I_pre(i-1)*s0*(1-K_pre(i));
    
    cum_I_pre(i)=cum_I_pre(i-1)+dI_pre(i);
    cum_R_pre(i)=cum_R_pre(i-1)+dR_pre(i);
    cum_D_pre(i)=cum_D_pre(i-1)+dD_pre(i);
    
    I_pre(i)=I_pre(i-1)+dI_pre(i)-dR_pre(i)-dD_pre(i);
end

C_pre = c0.*(1-G_pre);
S_pre = s0.*(1-K_pre);

yhat=log([I_pre;cum_I_pre;cum_R_pre;cum_D_pre;...
    dI_pre;dR_pre;dD_pre;...
    C_pre;Z_pre;S_pre]);
end
