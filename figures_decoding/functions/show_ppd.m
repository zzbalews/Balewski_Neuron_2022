function show_ppd(bhvdata,valdata,keep_tr, trialperiod)
%%
subject_names = unique(valdata.subject);
nsubj = length(subject_names);


% useful vars
t_mids = valdata.t_mids;


figure;
for s = 1:nsubj
   
    % restrict trials
    subject = subject_names{s};
    idx = find(strcmp(bhvdata.subject,subject) & ...
        keep_tr);
    
    rt = bhvdata.rt(idx);
    rt_cutoff = quantile(rt, 0.75);
    idx = idx(rt>=rt_cutoff);
    
    ntr = length(idx);
    
    % plot
    subplot(1,2,s); hold on
    
    if strcmp(trialperiod,'pics')
    rectangle('Position',[100 -20 400 200],'FaceColor',[.5 .5 .5 .2],'EdgeColor',[0 0 0 0])
    end
    
    ch = valdata.ch_ppd(idx,:);
    a = shadedErrorBar(t_mids, mean(ch),std(ch)./sqrt(ntr),...
        'lineprops',{'Color',[1 0 0]});
    
    unch = valdata.unch_ppd(idx,:);
    b = shadedErrorBar(t_mids, mean(unch), std(unch)./sqrt(ntr), ...
        'lineprops',{'Color',[0 0 1]});
    
    temp = randperm(ntr);
    ntr_half = round(ntr/2);
    na = cat(1,valdata.na_ppd(idx(temp(1:ntr_half)),:,1), ...
        valdata.na_ppd(idx(temp(ntr_half+1:end)),:,2));
    c = shadedErrorBar(t_mids, mean(na), std(na)./sqrt(ntr),...
        'lineprops',{'Color',[.5 .5 .5]});
    
%     na1 = valdata.na_ppd(idx,:,1);
%     shadedErrorBar(t_mids, mean(na1), std(na1)./sqrt(ntr),...
%         'lineprops',{'Color',[.5 .5 .5]})
%     na2 = valdata.na_ppd(idx,:,2);
%     shadedErrorBar(t_mids, mean(na2), std(na2)./sqrt(ntr),...
%         'lineprops',{'Color',[.5 .5 .5]})
%     plot(t_mids, mean(valdata.na_ppd(idx,:,1)),'Color',[.5 .5 .5])
%     plot(t_mids, mean(valdata.na_ppd(idx,:,2)),'Color',[.5 .5 .5])
    % format
    if strcmp(trialperiod,'pics')
    xlim([-500 1500])
    xlabel('time from pics (ms)')
    elseif strcmp(trialperiod,'choice')
        xlim([-1000 1000])
        xlabel('time from choice (ms)')
    end
    ylim([-10 150])
    
    set(gca,'XTick',-500:500:1500,'XTickLabel',-.5:.5:1.5,...
        'YTick',0:50:150)
    plot([-2000 2000],[0 0],'k:')
    plot([0 0],[-20 180],'k:')
    
    
    ylabel({'value decoder','%\Delta posterior probability'})
    title(subject)
    
    legend([a.mainLine,b.mainLine,c.mainLine],{'ch','unch','NA'},'box','off','Location','northwest')
    
    disp([subject,': peak value decoding'])
    [~,x_ch] = max(mean(ch(:,t_mids<500)));
    disp(['   ch: ',num2str(1+t_mids(x_ch)),'ms'])
    [~,x_unch] = max(mean(unch(:,t_mids<500)));
    disp(['   unch: ',num2str(1+t_mids(x_unch)),'ms'])
end


end