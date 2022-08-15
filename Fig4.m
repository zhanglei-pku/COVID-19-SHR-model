% This file is used to plot the figure 4
% Author: Lei Zhang
% Last modified: 2022-06-27


clearvars;
close all;
clc;

% the bound of parameters of separating degree (i0, gamma, ts, s1)
S_lb = [0    0   0    0];
S_ub = [1    1   1000  1];

% the bound of parameters of healing degree (alpha, h1, th)
H_lb = [0    0   0 ];
H_ub = [1    1   1000];

tableSD = readtable('province.csv');
allprovince = {...
    'Anhui';...
    'Beijing';...
    'Chongqing';...
    'Fujian';...
    'Guangdong';...
    'Guangxi';...
    'Guizhou';...
    'Hainan';...
    'Hebei';...
    'Heilongjiang';...
    'Henan';...
    'Hunan';...
    'Jiangsu';...
    'Jiangxi';...
    'Liaoning';...
    'Shaanxi';...
    'Shandong';...
    'Shanghai';...
    'Shanxi';...
    'Sichuan';...
    'Tianjin';...
    'Yunnan';...
    'Zhejiang'};

gcf = figure('position',[284,121,1065,755]);
gcf_t = tiledlayout(5,5,'TileSpacing','tight','Padding','tight');
for k = 1 :length(allprovince)
    clearvars -except k allprovince tableSD S_lb S_ub H_lb H_ub gcf gcf_t
    
    inputarea = allprovince{k};
    if isempty(inputarea)
       warning('Please enter the country or region and then press the confirmation button ！')
       return
    end

    %% Download the data from ref [1] and read them with the function func_getDataCOVID

    [tableConfirmed,tableDeaths,tableRecovered,time] = func_getDataCOVID();
    fprintf(['Most recent update: ',datestr(time(end)),'\n'])

    if ~isempty(find(tableRecovered.ProvinceState==inputarea, 1))
        indI = find(tableConfirmed.ProvinceState==inputarea);
        indC = find(tableRecovered.ProvinceState==inputarea);   
        indD = find(tableDeaths.ProvinceState==inputarea);
        
        Confirmed = table2array(tableConfirmed(indI,5:end));
        Recovered = table2array(tableRecovered(indC,5:end));
        Deaths    = table2array(tableDeaths(indD,5:end));
        
        disp(tableConfirmed(indI(1),1:2))
        disp(tableRecovered(indC(1),1:2))
        disp(tableDeaths(indD(1),1:2))
    else
        % Discuss the different situations of the Confirmed

        if isempty(find(tableConfirmed.CountryRegion==inputarea, 1))
            warning('Could not find the country or region, please check the inputarea (The first letter should be capitalized).')
            return
        elseif ~isempty(find((tableConfirmed.CountryRegion==inputarea) & (tableConfirmed.ProvinceState.ismissing()==1), 1))
            indI = find((tableConfirmed.CountryRegion==inputarea) & (tableConfirmed.ProvinceState.ismissing()==1));
            Confirmed = table2array(tableConfirmed(indI,5:end));
            disp(tableConfirmed(indI(1),1:2))
        else
            indI = find((tableConfirmed.CountryRegion==inputarea));
            Confirmed = sum(table2array(tableConfirmed(indI,5:end)),1);
            disp(tableConfirmed(indI(1),2:2))
        end
        
        % Discuss the different situations of the Recovered
        if isempty(find(tableRecovered.ProvinceState==inputarea, 1))
            warning('Could not find the country or region, please check the inputarea. The first letter should be capitalized.')
            return
        elseif ~isempty(find((tableRecovered.ProvinceState==inputarea) & (tableRecovered.ProvinceState.ismissing()==1), 1))
            indC = find((tableRecovered.ProvinceState==inputarea) & (tableRecovered.ProvinceState.ismissing()==1));
            Recovered = table2array(tableRecovered(indC,5:end));
            disp(tableRecovered(indC(1),1:2))
        else
            indC = find((tableRecovered.ProvinceState==inputarea));
            Recovered = sum(table2array(tableRecovered(indC,5:end)),1);
            disp(tableRecovered(indC(1),2:2))
        end

        % Discuss the different situations of the Deaths

        if isempty(find(tableDeaths.CountryRegion==inputarea, 1))
            warning('Could not find the country or region, please check the inputarea (The first letter should be capitalized).')
            return
        elseif ~isempty(find((tableDeaths.CountryRegion==inputarea) & (tableDeaths.ProvinceState.ismissing()==1), 1))
            indD = find((tableDeaths.CountryRegion==inputarea) & (tableDeaths.ProvinceState.ismissing()==1));
            Deaths = table2array(tableDeaths(indD,5:end));
            disp(tableDeaths(indD(1),1:2))
        else
            indD = find((tableDeaths.CountryRegion==inputarea));
            Deaths = sum(table2array(tableDeaths(indD,5:end)),1);
            disp(tableDeaths(indD(1),2:2))
        end
    end

    %% Prepare the data from the date when having the cases 

    ind_zero_Confirmed = find(Confirmed <= 0);

    if ~isempty(ind_zero_Confirmed)
        Confirmed = Confirmed(ind_zero_Confirmed(end)+1:end);
        Recovered = Recovered(ind_zero_Confirmed(end)+1:end);
        Deaths    = Deaths(ind_zero_Confirmed(end)+1:end);
        time = time(ind_zero_Confirmed(end)+1:end);
    end

    if isempty(Confirmed)
        warning('Confirmed is empty.')
        return
    end
    
    
    %% prepare the data to be fitted
    % the startpoint and endpoint
    indSD   = find(strcmp(tableSD.location,inputarea)==1);
    
    ind_start_confirmed = tableSD.ind_start1_confirmed(indSD);
    ind_end_confirmed   = tableSD.ind_end1_confirmed(indSD);
    
    ind_start_recovered = tableSD.ind_start1_recovered(indSD);
    ind_end_recovered   = tableSD.ind_end1_recovered(indSD);
    
    % the time of prediction
    t = 1:length(time);
    time_pre = time(1:200);
    t_confirmed = t(ind_start_confirmed:ind_end_confirmed);
    t_recovered = t(ind_start_recovered:ind_end_recovered);
    t_quarantined = t(ind_start_confirmed:ind_end_recovered);
    
    % cumulative confirmed
    cum_confirmed_fit = Confirmed(ind_start_confirmed:ind_end_confirmed);
    
    % new confirmed
    new_confirmed0 = [0,diff(cum_confirmed_fit)];
    index_newconf = ~isinf(log(new_confirmed0));
    t_newconf = t_confirmed(index_newconf);
    t_newconf_spl = t_newconf(1):1:t_newconf(end);
    new_confirmed_fit = interp1(t_newconf,new_confirmed0(index_newconf),...
        t_newconf_spl,'makima');
    new_confirmed_fit(new_confirmed_fit<0) = 1;
    
    % cumulative recovered
    cum_recovered_fit = Recovered(ind_start_recovered:ind_end_recovered);
    
    % new recovered
    new_recovered0 = [0,diff(cum_recovered_fit)];
    index_newreco = ~isinf(log(new_recovered0));
    t_newreco = t_recovered(index_newreco);
    t_newreco_spl = t_newreco(1):1:t_newreco(end);
    new_recovered_fit = interp1(t_newreco,new_recovered0(index_newreco),...
        t_newreco_spl,'makima');
    new_recovered_fit(new_recovered_fit<0) = 1;
    
    % quarantined
    quarantined_fit = Confirmed(ind_start_confirmed:ind_end_recovered) - ...
        Recovered(ind_start_confirmed:ind_end_recovered);

    % date
    date_cum_confirmed = time(t_confirmed);
    date_cum_recovered = time(t_recovered);

    date_new_confirmed = time(t_newconf_spl);
    date_new_recovered = time(t_newreco_spl);

    date_quarantined = time(t_quarantined);
    
    % the fitted data
    DataTobeFitted=transpose([...
        quarantined_fit,...
        cum_confirmed_fit,...
        cum_recovered_fit,...
        new_confirmed_fit,...
        new_recovered_fit]);
    
    % the length of fitted data
    N1=length(quarantined_fit);
    N2=length(cum_confirmed_fit);
    N3=length(cum_recovered_fit);
    N4=length(new_confirmed_fit);
    N5=length(new_recovered_fit);
    N_all =[N1,N2,N3,N4,N5];
    
    % the transition rate
    new_confirmed_fit2 = [0,new_confirmed_fit];
    new_recovered_fit2 = [0,new_recovered_fit];

    for i = 2:N_all(4)+1
        rate_qz(i)   = new_confirmed_fit2(i)./quarantined_fit(i-1);
    end

    for i = 2:N_all(5)+1
        rate_cure(i) = new_recovered_fit2(i)./quarantined_fit(i-1);
    end

    %% the first fitting: fit the transtion rate
    
    SetPara0=[0.5 0.5 0.5 0.5];
    [paraC,ypreC1] = func_fit_S(log(rate_qz(2:end)),SetPara0,...
        t_newconf_spl,S_lb,S_ub);

    SetPara0=[0.1 0.5 1];
    [paraZ,ypreZ1] = func_fit_H(log(rate_cure(2:end)),SetPara0,...
        t_newreco_spl,H_lb,H_ub) ;

    paraRate = [paraC,paraZ]';
    Inputpara0=[DataTobeFitted(1);paraRate];
    
    t_gap_cumcure = t_recovered(1);
    t_gap_newcure = t_newreco_spl;
    
    [RSD0,ypre0] = func_simulate_rate(...
        DataTobeFitted,Inputpara0(:,1),N_all,t_gap_cumcure,t_gap_newcure);

    %% the second fitting: fit the cases data
    
    
    % Call the fitting program
    ErrorCriterion=10^(-24);
    [paraIRD_Optimal,RSD_Optimal,ypre_Optimal] = func_fit_cases(...
        DataTobeFitted,N_all,RSD0,Inputpara0,ErrorCriterion,...
        t_gap_cumcure,t_gap_newcure,S_lb,S_ub,H_lb,H_ub);

    num_length =length(DataTobeFitted)+N_all(4)+N_all(5);

    for i =1:9
        data_preout(:,i)= exp(ypre_Optimal(num_length+(i-1)*200+1:num_length+i*200));
    end

    now_qz_preout = data_preout(:,1);
    cum_qz_preout = data_preout(:,2);
    cum_cure_preout = data_preout(:,3);
    new_qz_preout = data_preout(:,4);
    new_cure_preout = data_preout(:,5);

    I0_preout = paraIRD_Optimal(1);
    i0_preout = paraIRD_Optimal(2);
    gamma_preout = paraIRD_Optimal(3);
    ts_preout = paraIRD_Optimal(4);
    s1_preout = paraIRD_Optimal(5);

    alpha_preout = paraIRD_Optimal(6);
    h1_preout = paraIRD_Optimal(7);
    th_preout = paraIRD_Optimal(8);
    
    S1_preout = s1_preout/(1+exp(-gamma_preout*(1-ts_preout)));
    H1_preout = h1_preout/(1+exp(-alpha_preout*(1-th_preout)));

    %% duration
    [new_qz_max,t_qz_max] = max(new_qz_preout);
    ind_qz_1   = find(new_qz_preout(2:t_qz_max)-new_qz_max/10>0);
    ind_qz_2   = find(new_qz_preout(t_qz_max:end)-new_qz_max/10<0);
    if isempty(ind_qz_2)
        duration_qz  = 0;
    else
        t_qz_1 = ind_qz_1(1)+1;
        t_qz_2 = ind_qz_2(1)-1+t_qz_max - 1;
        duration_qz = t_qz_2 - t_qz_1 + 1;
    end


    [new_cure_max,t_cure_max] = max(new_cure_preout);
    ind_cure_1   = find(new_cure_preout(2:t_cure_max)-new_cure_max/10>0);
    ind_cure_2   = find(new_cure_preout(t_cure_max:end)-new_cure_max/10<0);
    if isempty(ind_cure_2)
        duration_cure  = 0;
    else
        t_cure_1 = ind_cure_1(1)+1;
        t_cure_2 = ind_cure_2(1)-1+t_cure_max - 1;
        duration_cure   = t_cure_2 - t_cure_1 + 1;
    end


    [now_qz_max,t_now_max] = max(now_qz_preout);
    ind_now_1   = find(now_qz_preout(2:t_now_max)-now_qz_max/10>0);
    ind_now_2   = find(now_qz_preout(t_now_max:end)-now_qz_max/10<0);
    if isempty(ind_now_2)
        duration_now   = 0;
    else
        t_now_1 = ind_now_1(1)+1;
        t_now_2 = ind_now_2(1)-1+t_now_max - 1;
        duration_now   = t_now_2 - t_now_1 + 1;
    end
    
    
    %% plot the results
    
    func_plot_all_province(DataTobeFitted,data_preout,N_all,time_pre,...、
        date_cum_confirmed,date_cum_recovered,...
        date_new_confirmed,date_new_recovered,...
        date_quarantined,inputarea);
    
    
    
    % caculte the key date
    date_1= datenum(time(1));

    % ts
    date_ts = datestr(date_1+ts_preout,'yyyy-mm-dd');
    date_ts = datetime(date_ts,'InputFormat','yyyy-MM-dd');
    ts_new = date_1+ts_preout - datenum('2020-01-20');

    % th
    date_th = datestr(date_1+th_preout,'yyyy-mm-dd');
    date_th =datetime(date_th,'InputFormat','yyyy-MM-dd');
    th_new = date_1+th_preout - datenum('2020-01-20');
    
    % save the results
    tableSD.ind_start_confirmed(indSD) = ind_start_confirmed;
    tableSD.ind_end_confirmed(indSD)   = ind_end_confirmed;
    tableSD.ind_start_recovered(indSD) = ind_start_recovered;
    tableSD.ind_end_recovered(indSD)   = ind_end_recovered;
    
    tableSD.c0(indSD) = i0_preout;
    tableSD.gamma(indSD) = gamma_preout;
    tableSD.g0(indSD) = s1_preout;
    tableSD.t_g(indSD) = ts_preout;
    
    tableSD.alpha(indSD) = alpha_preout;
    tableSD.h0(indSD)  = h1_preout;
    tableSD.t_h(indSD) = th_preout;
    
    
    tableSD.date_tg(indSD)= date_ts;
    tableSD.tg_new(indSD) = ts_new;
    tableSD.date_th(indSD)= date_th;
    tableSD.th_new(indSD) = th_new;
    
    tableSD.cum_cases(indSD) = Confirmed(ind_end_confirmed);
    tableSD.RSD(indSD) = RSD_Optimal;
    
    tableSD.t_qz(indSD)   = duration_qz;
    tableSD.t_cure(indSD) = duration_cure;
    tableSD.t_now(indSD)  = duration_now;
    
    tableSD.G1_pre(indSD)  = S1_preout;
    tableSD.H1_pre(indSD)  = H1_preout;
  
    
end


lgd = legend({'daily infected','daily cured','total infected','total cured',...
    'model'},'FontSize',12,'numcolumns',2);
lgd.Layout.Tile = [24,25];
lad.Position =[0.672707356771393,0.068344370083304,0.254460089167519,0.090066222678747];
% title(gcf_t,'23 areas in the first wave','FontWeight','bold','FontSize',16)

% writetable(tableSD,'provence.csv');





%% user-defined functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fit the infection rate
function [para,ypre] = func_fit_S(rate,SetPara0,t_newconf_spl,S_ul,S_ub)
tData = t_newconf_spl;
yData = rate;

lb = S_ul;
ub = S_ub;

% options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
a = lsqcurvefit(@AnalyticFun,SetPara0,tData,yData,lb,ub');
ypre=AnalyticFun(a,tData);
para=a;

function yhat = AnalyticFun(a,t)
Vi0 = a(1);
gamma = a(2);
ts = a(3);
s1 = a(4);

St = s1./(1+exp(-gamma*(t-ts)));

Vi = Vi0*(1-St);

yhat = log(Vi);

end

end


% fit the cure rate
function [para,ypre] = func_fit_H(rate,SetPara0,t_newreco_spl,H_ul,H_ub)
tData = t_newreco_spl;
yData = rate;

lb = H_ul;
ub = H_ub;
a = lsqcurvefit(@AnalyticFun,SetPara0,tData,yData,lb,ub);
ypre=AnalyticFun(a,tData);
para=a;

function yhat=AnalyticFun(a,t)
alpha = a(1);
h1 = a(2);
th = a(3);

Vc = h1./(1+exp(-alpha*(t-th)));

yhat=log(Vc);

end

end

% simulate the transition rate
function [RSD,ypre] = func_simulate_rate(Data,Inputpara,N_all,t_gap_cumcure,t_gap_newcure)

cum_I1 = Data(N_all(1)+1,1);
cum_C1 = Data(N_all(1)+N_all(2)+1,1);

new_I1 = Data(N_all(1)+N_all(2)+N_all(3)+1,1);
new_C1 = Data(N_all(1)+N_all(2)+N_all(3)+N_all(4)+1,1);

yData = log(Data);

ypre = OdeNew(N_all,Inputpara,cum_I1,cum_C1,new_I1,new_C1,t_gap_cumcure,t_gap_newcure);
residual = ypre(1:length(yData))-yData;
SSE = sum(residual.^2);
RSD = sqrt(SSE/(length(yData)-2));

end

% modle
function yhat = OdeNew(N_all,Inputpara,cum_I1,cum_C1,new_I1,new_C1,t_gap_cumcure,t_gap_newcure)


I0 = Inputpara(1);
Vi0 = Inputpara(2);
gamma = Inputpara(3);
ts = Inputpara(4);
s1 = Inputpara(5);

alpha = Inputpara(6); 
h1 = Inputpara(7);      
th = Inputpara(8);


N_max = max(N_all);

I_pre = zeros(N_max,1);
S_pre = zeros(N_max,1);
H_pre = zeros(N_max,1);

dI_pre = zeros(N_max,1);
dC_pre = zeros(N_max,1);

cum_I_pre = zeros(N_max,1);
cum_C_pre = zeros(N_max,1);

S0 = s1/(1+exp(-gamma*(1-ts)));
H0 = h1/(1+exp(-alpha*(1-th)));

I_pre(1) = I0;
cum_I_pre(1) = cum_I1;

S_pre(1) = S0;
H_pre(1) = H0;

dI_pre(1) = 0;
dI_pre(2) = new_I1;



for i = 1:t_gap_cumcure
    cum_C_pre(i) = 0;
end


for i = 2:1:200
    S_pre(i) = s1/(1+exp(-gamma*(i-ts)));
    H_pre(i) = h1/(1+exp(-alpha*(i-th)));
    
    dI_pre(i)= I_pre(i-1)*Vi0*(1-S_pre(i));
    dC_pre(i)= I_pre(i-1)*H_pre(i);

    cum_I_pre(i) = cum_I_pre(i-1) + dI_pre(i);
    cum_C_pre(i) = cum_C_pre(i-1) + dC_pre(i);

    I_pre(i) = I_pre(i-1) + dI_pre(i) - dC_pre(i);
end
Vi_pre = Vi0.*(1-S_pre);
Vc_pre = H_pre;

yhat=log([...
    I_pre(1:N_all(1));...
    cum_I_pre(1:N_all(2));...
    cum_C_pre(t_gap_cumcure:(N_all(3)+t_gap_cumcure-1));...
    dI_pre(2:N_all(4)+1);...
    dC_pre(t_gap_newcure:(N_all(5)+t_gap_newcure-1));...
    Vi_pre(1:N_all(4));...
    Vc_pre(1:N_all(5));...
    I_pre;...
    cum_I_pre;...
    cum_C_pre;...
    dI_pre;...
    dC_pre;...
    Vi_pre;...
    Vc_pre;...
    S_pre;...
    H_pre;]);

end
