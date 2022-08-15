function [DataTobeFitted,data_pre2,N_all,time_pre,...
    date_cum_confirmed,date_cum_recovered,...
    date_new_confirmed,date_new_recovered,...
    date_quarantined,paraIRD_Optimal]  = func_fit_3wave(kk,inputarea)

% the parameters of separating degree S(i0, gamma, ts, s1)
S_lb = [0     0    0     0];
S_ub = [3     1    100   1];

% the parameters of healing degree H(alpha, h1, th)
H_lb = [0   0    0];
H_ub = [1   0.5    100];

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
    
    Recovered = table2array(tableRecovered(indC,5:end));
    Deaths    = table2array(tableDeaths(indD,5:end));
    Confirmed = table2array(tableConfirmed(indI,5:end));
    
    disp(tableRecovered(indC(1),1:2))
    disp(tableConfirmed(indI(1),1:2))
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
    Recovered = Recovered(ind_zero_Confirmed(end)+1:end);
    Deaths    = Deaths(ind_zero_Confirmed(end)+1:end);
    Confirmed = Confirmed(ind_zero_Confirmed(end)+1:end);
    time = time(ind_zero_Confirmed(end)+1:end);
end

if isempty(Confirmed)
    warning('Confirmed is empty.')
    return
end


t = 1:length(time);



%% Determine the endpoint of Confirmed

tableSD = readtable('province.csv');
indSD   = find(strcmp(tableSD.location,inputarea)==1);



ind_start3_confirmed = tableSD.ind_start3_confirmed(indSD);
ind_end3_confirmed   = tableSD.ind_end3_confirmed(indSD)-kk+1;

ind_start3_recovered  = tableSD.ind_start3_recovered(indSD);
ind_end3_recovered  =   tableSD.ind_end3_recovered(indSD)-kk+1;


%% Caculate the daily cases

t_confirmed = t(ind_start3_confirmed:ind_end3_confirmed);
t_recovered = t(ind_start3_recovered:ind_end3_recovered);
t_quarantined = t(ind_start3_confirmed:ind_end3_recovered);


% confirmed
cum_confirmed_fit = Confirmed(ind_start3_confirmed:ind_end3_confirmed);
cum_confirmed_fit = cum_confirmed_fit - Confirmed(ind_start3_confirmed-1);

new_confirmed0 = [0,diff(cum_confirmed_fit)];


% recovered
cum_recovered_fit = Recovered(ind_start3_recovered:ind_end3_recovered);
cum_recovered_fit = cum_recovered_fit - Recovered(ind_start3_recovered-1);
new_recovered0 = [0,diff(cum_recovered_fit)];

% Quarantined are the cases that are still treated in hospital
gap_data = zeros(1,length(cum_confirmed_fit)-length(cum_recovered_fit));
quarantined_fit = cum_confirmed_fit - [gap_data,cum_recovered_fit];

% new confirmed
index_newconf = ~isinf(log(new_confirmed0));
t_newconf = t_confirmed(index_newconf);
t_newconf_spl = t_newconf(1):1:t_newconf(end);
new_confirmed_fit = interp1(t_newconf,new_confirmed0(index_newconf),...
    t_newconf_spl,'makima');
new_confirmed_fit(new_confirmed_fit<0) = 1;

% new recovered
index_newreco = ~isinf(log(new_recovered0));
t_newreco = t_recovered(index_newreco);
t_newreco_spl = t_newreco(1):1:t_newreco(end);
new_recovered_fit = interp1(t_newreco,new_recovered0(index_newreco),...
    t_newreco_spl,'makima');%
new_recovered_fit(new_recovered_fit<0) = 1;


%% date
date_cum_confirmed = time(t_confirmed);
date_cum_recovered = time(t_recovered);

date_new_confirmed = time(t_newconf_spl);
date_new_recovered = time(t_newreco_spl);

date_quarantined = time(t_quarantined);

time_pre =  time(t_confirmed(1)):1:time(t_confirmed(1))+199;


%% remove the 0 or NaN data

N1=length(quarantined_fit);
N2=length(cum_confirmed_fit);
N3=length(cum_recovered_fit);
N4=length(new_confirmed_fit);
N5=length(new_recovered_fit);
N_all =[N1,N2,N3,N4,N5];

DataTobeFitted=transpose([...
    quarantined_fit,...
    cum_confirmed_fit,...
    cum_recovered_fit,...
    new_confirmed_fit,...
    new_recovered_fit]);

%% rate
new_confirmed_fit2 = [0,new_confirmed_fit];
new_recovered_fit2 = [0,new_recovered_fit];

for i = 2:N_all(4)+1
    rate_qz(i)   = new_confirmed_fit2(i)./quarantined_fit(i-1);
end

for i = 2:N_all(5)+1
    rate_cure(i) = new_recovered_fit2(i)./quarantined_fit(i-1);
end


%% the first step of fitting: fit the rate data.
% infection rate
SetPara0=[0.5 0.5 0.5 0.5];
[paraC,ypreC1]= nlinfit_C_Simple(log(rate_qz(2:end)),SetPara0,...
    t_newconf_spl,S_lb,S_ub);

% cure rate
SetPara0=[0.1 0.5 1];
[paraZ,ypreZ1] =nlinfit_Z_Simple(log(rate_cure(2:end)),SetPara0,...
    t_newreco_spl,H_lb,H_ub) ;

paraRate = [paraC,paraZ]';
Inputpara0=[DataTobeFitted(1);paraRate];

t_gap_cumcure = t_recovered(1)-t_confirmed(1)+1;
t_gap_newcure = t_newreco_spl(1)-t_confirmed(1)+1;

[RSD0,ypre0] = fit_result_model(DataTobeFitted,Inputpara0(:,1),N_all,t_gap_cumcure,t_gap_newcure);



%% the second step of fitting: fit the case data.


ErrorCriterion=10^(-6);
[paraIRD_Optimal,RSD_Optimal,ypre_Optimal]=func_fit_cases(DataTobeFitted,N_all,...
    RSD0,Inputpara0,ErrorCriterion,t_gap_cumcure,t_gap_newcure,S_lb,S_ub,H_lb,H_ub);

num_length =length(DataTobeFitted)+N_all(4)+N_all(5);

for i =1:9
    data_pre2(:,i)= exp(ypre_Optimal(num_length+(i-1)*200+1:num_length+i*200));
end



%% user-defined functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% fit the infection rate
function [para,ypre]=nlinfit_C_Simple(rate,SetPara0,t_newconf_spl,S_lb,S_ub)
tData = t_newconf_spl;
yData = rate;

lb = S_lb;
ub = S_ub;
% options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
a = lsqcurvefit(@AnalyticFun,SetPara0,tData,yData,lb,ub');
ypre=AnalyticFun(a,tData);
para=a;
% figure
% hold on
% plot(tData,yData,'ro')
% plot(tData,ypre,'b_')
% plot(tData,ypre,'b_')

function yhat=AnalyticFun(a,t)
Vi0 = a(1);
gamma = a(2);
ts = a(3);
s1 = a(4);

St = s1./(1+exp(-gamma*(t-ts)));

Vi = Vi0*(1-St);

yhat=log(Vi);

end

end


% fit the cure rate
function [para,ypre] = nlinfit_Z_Simple(rate,SetPara0,t_newreco_spl,H_lb,H_ub)
tData = t_newreco_spl;
yData = rate;

tData2=-100:1:100;

lb = H_lb;
ub = H_ub;
a = lsqcurvefit(@AnalyticFun,SetPara0,tData,yData,lb,ub);
ypre=AnalyticFun(a,tData2);
para=a;

function yhat=AnalyticFun(a,t)
alpha = a(1);
h1 = a(2);
th = a(3);

Vc = h1./(1+exp(-alpha*(t-th)));

yhat=log(Vc);

end

end

% calculate the cases
function [RSD,ypre] = fit_result_model(Data,Inputpara,N_all,t_gap_cumcure,t_gap_newcure)

I0     = Data(1);
cum_I1 = Data(N_all(1)+1,1);
cum_C1 = Data(N_all(1)+N_all(2)+1,1);

new_I1 = Data(N_all(1)+N_all(2)+N_all(3)+1,1);
new_C1 = Data(N_all(1)+N_all(2)+N_all(3)+N_all(4)+1,1);

yData = log(Data);

ypre = OdeNew(N_all,Inputpara,cum_I1,new_I1,t_gap_cumcure,t_gap_newcure);
residual = ypre(1:length(yData))-yData;
SSE = sum(residual.^2);
RSD = sqrt(SSE/(length(yData)-2));

end

function yhat = OdeNew(N_all,Inputpara,cum_I1,new_I1,t_gap_cumcure,t_gap_newcure)


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


cum_I_pre(1) = cum_I1;
I_pre(1) = I0;
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
end
