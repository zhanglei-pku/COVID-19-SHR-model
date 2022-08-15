

function [Outputpara,RSD,ypre_best]=func_fit_cases(Data,N_all,...
    RSD0,Inputpara0,ErrorCriterion,t_gap_cumcure,t_gap_newcure,S_lb,S_ub,H_lb,H_ub)

    tData=1:1:length(Data);

    cum_I1=Data(N_all(1)+1,1);
    cum_C1=Data(N_all(1)+N_all(2)+1,1);
    new_I1=Data(N_all(1)+N_all(2)+N_all(3)+1,1);
    new_C1=Data(N_all(1)+N_all(2)+N_all(3)++N_all(4)+1,1);

    yData=log(Data);

    
    flag_xunhuan=1;
    num_error =1;
    min_error(1) = RSD0;
    InputparaZ=Inputpara0;
    while flag_xunhuan==1
        num_error = num_error+1;
        [InputparaC,RSD1,ypre1]=fitC(tData,yData,InputparaZ,t_gap_cumcure,...
            t_gap_newcure,S_lb,S_ub);
        
        [InputparaZ,RSD2,ypre2]=fitZ(tData,yData,InputparaC,t_gap_cumcure,...
            t_gap_newcure,H_lb,H_ub);
        
        if RSD2 > RSD1
            InputparaZ = InputparaC;
            RSD = RSD1;
            ypre_best = ypre1;
        end
        
        RSD = RSD2;
        ypre_best = ypre2;
        min_error(num_error)=RSD;
        
        if abs(min_error(num_error-1)-min_error(num_error))<ErrorCriterion...
                || num_error > 100
            flag_xunhuan = 2;
        else
            flag_xunhuan = 1;
        end
    end
    Outputpara = InputparaZ;

function [Inputpara,RSD,ypre1]=fitC(tData,yData,Inputpara,t_gap_cumcure,...
        t_gap_newcure,S_ul,S_ub)

    SetPara0=Inputpara0(2:5);

    lb = S_ul;
    ub = S_ub;
    a = lsqcurvefit(@CFun,SetPara0,tData,yData,lb,ub);

    Inputpara(2)=a(1);
    Inputpara(3)=a(2);
    Inputpara(4)=a(3);
    Inputpara(5)=a(4);

    ypre = CFun(a,tData);
    residual = ypre(1:length(yData))-yData;
    SSE = sum(residual.^2);
    RSD = sqrt(SSE/(length(yData)-2)); 
    ypre1 = CFun2(a,tData);

    function yhat2=CFun(a,tData)
        
        I0 = Inputpara(1);
        Vi0 = a(1);
        gamma = a(2);
        ts = a(3);
        s1 = a(4);

        alpha = Inputpara(6); 
        h1 = Inputpara(7);    
        th = Inputpara(8);

        yhat1 = OdeNew(N_all, Vi0,gamma,ts,s1,alpha,h1,th,I0,...
            cum_I1,cum_C1,new_I1,new_C1,t_gap_cumcure,t_gap_newcure);
        yhat2 = yhat1(1:length(tData));
    end
    
    function yhat1=CFun2(a,tData)
        
        I0 = Inputpara(1);
        Vi0 = a(1);
        gamma = a(2);
        ts = a(3);
        s1 = a(4);

        alpha = Inputpara(6); 
        h1 = Inputpara(7);    
        th = Inputpara(8);

        yhat1 = OdeNew(N_all, Vi0,gamma,ts,s1,alpha,h1,th,I0,...
            cum_I1,cum_C1,new_I1,new_C1,t_gap_cumcure,t_gap_newcure);
    end

end

    function [Inputpara,RSD,ypre1]=fitZ(tData,yData,Inputpara,...
            t_gap_cumcure,t_gap_newcure,H_ul,H_ub)
    SetPara0=Inputpara0(6:8);

    lb = H_ul;
    ub = H_ub;
    a = lsqcurvefit(@ZFun,SetPara0,tData,yData,lb,ub);
    Inputpara(6) = a(1);
    Inputpara(7) = a(2);
    Inputpara(8) = a(3);

    ypre = ZFun(a,tData);
    residual = ypre(1:length(yData))-yData;
    SSE = sum(residual.^2);
    RSD = sqrt(SSE/(length(yData)-2)); 
    ypre1 = ZFun2(a,tData);

    function yhat2=ZFun(a,tData)
        
        I0 = Inputpara(1);
        Vi0 = Inputpara(2);
        gamma = Inputpara(3);
        ts = Inputpara(4);
        s1 = Inputpara(5);

        alpha = a(1); 
        h1 = a(2);    
        th = a(3);

        yhat1 = OdeNew(N_all,Vi0,gamma,ts,s1,alpha,h1,th,I0,...
            cum_I1,cum_C1,new_I1,new_C1,t_gap_cumcure,t_gap_newcure);
        yhat2 = yhat1(1:length(tData));
    end
    
    function yhat1 = ZFun2(a,tData)
        
        I0 = Inputpara(1);
        Vi0 = Inputpara(2);
        gamma = Inputpara(3);
        ts = Inputpara(4);
        s1 = Inputpara(5);

        alpha = a(1); 
        h1 = a(2);    
        th = a(3);

        yhat1=OdeNew(N_all,Vi0,gamma,ts,s1,alpha,h1,th,I0,...
            cum_I1,cum_C1,new_I1,new_C1,t_gap_cumcure,t_gap_newcure);
    end
end



function yhat=OdeNew(N_all,Vi0,gamma,ts,s1,alpha,h1,th,...
        I0,cum_I1,cum_C1,new_I1,new_C1,t_gap_cumcure,t_gap_newcure)
    
    N_max = max(N_all);

    I_pre=zeros(N_max,1);
    S_pre=zeros(N_max,1);
    H_pre=zeros(N_max,1);
    dI_pre=zeros(N_max,1);
    dC_pre=zeros(N_max,1);
    cum_I_pre=zeros(N_max,1);
    cum_C_pre=zeros(N_max,1);

    S0 = s1/(1+exp(-gamma*(1-ts)));
    H0 = h1/(1+exp(-alpha*(1-th)));
    
    cum_I_pre(1) = cum_I1;
    I_pre(1) = I0;
    S_pre(1) = S0;
    H_pre(1) = H0;

    for i = 1:t_gap_cumcure
        cum_C_pre(i) = 0;
    end

    dI_pre(1) = 0;
    dI_pre(2) = new_I1;


    for i=2:1:200
        S_pre(i) = s1/(1+exp(-gamma*(i-ts)));
        H_pre(i) = h1/(1+exp(-alpha*(i-th)));

        dI_pre(i)=I_pre(i-1)*Vi0*(1-S_pre(i));
        dC_pre(i)=I_pre(i-1)*H_pre(i);

        cum_I_pre(i)=cum_I_pre(i-1)+dI_pre(i);
        cum_C_pre(i)=cum_C_pre(i-1)+dC_pre(i);

        I_pre(i)=I_pre(i-1)+dI_pre(i)-dC_pre(i);
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