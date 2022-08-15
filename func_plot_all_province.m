function func_plot_all_province(DataTobeFitted,data_pre,N_all,time_pre,...
    date_cum_confirmed,date_cum_recovered,...
    date_new_confirmed,date_new_recovered,...
    date_quarantined,inputarea)

    

    quarantined_fit   = DataTobeFitted(1:N_all(1));
    cum_confirmed_fit = DataTobeFitted(N_all(1)+1:N_all(1)+N_all(2));
    cum_recovered_fit = DataTobeFitted(N_all(1)+N_all(2)+1:N_all(1)+N_all(2)+N_all(3));
    new_confirmed_fit = DataTobeFitted(N_all(1)+N_all(2)+N_all(3)+1:N_all(1)+N_all(2)+N_all(3)+N_all(4));
    new_recovered_fit = DataTobeFitted(N_all(1)+N_all(2)+N_all(3)+N_all(4)+1:N_all(1)+N_all(2)+N_all(3)+N_all(4)+N_all(5));

    
    cum_qz_pre   = data_pre(:,2);
    cum_cure_pre = data_pre(:,3);
    new_qz_pre   = data_pre(:,4);
    new_cure_pre = data_pre(:,5);

%% Plot
%     gcf = figure;
%     t = tiledlayout(1,1,'TileSpacing','loose','Padding','compact'); 

    nexttile([1,1])
    hold on
    
    
    % new case
    plot(date_new_confirmed,new_confirmed_fit,'^','MarkerSize',3,'linewidth',1,...
        'color',[0.80392	0.33333	0.33333]);
    plot(date_new_recovered,new_recovered_fit,'k^','MarkerSize',3,'linewidth',1,...
        'color',[0.27451	0.5098	0.70588]);
    
    % cumulative case
    plot(date_cum_confirmed,cum_confirmed_fit,'ro','MarkerSize',3,'linewidth',1,...
        'color',[0.80392	0.33333	0.33333]);
    plot(date_cum_recovered,cum_recovered_fit,'ko','MarkerSize',3,'linewidth',1,...
        'color',[0.27451	0.5098	0.70588]);
    
    % model
    plot(time_pre,cum_qz_pre,'k-','LineWidth',1)
    plot(time_pre,cum_cure_pre,'k-','LineWidth',1)
    plot(time_pre,new_qz_pre,'k-','LineWidth',1)
    plot(time_pre,new_cure_pre,'k-','LineWidth',1)
    
    
    ylabel('Cases','FontSize',12)
    set(gca,'yscale','log')
    set(gca,'FontSize',12)
    box on
    xlim([time_pre(1) time_pre(length(quarantined_fit)+10)])
    ylim([1,10^4])
    
    date_xtick = [date_cum_confirmed(1),...
        time_pre(length(quarantined_fit)+10)];
    set(gca,'xtick',date_xtick);
    
    
    title(inputarea,'FontSize',12)
    
end