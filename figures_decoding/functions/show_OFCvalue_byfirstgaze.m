function show_OFCvalue_byfirstgaze(bhvdata,decdata,tr_sampled,decoder,keep_tr)


%% useful vars
subject_names = unique(decdata.subject);
nsubj = length(subject_names);

t_mids = decdata.t_mids;
nmids = length(t_mids);
[~,t0] = min(abs(decdata.t_mids-200));
[~,t_end] = min(abs(decdata.t_mids-500));

niters = size(tr_sampled.(subject_names{1}).keep_match,2);

ntr = length(bhvdata.TrialNumber);
%% get first OFC value state after pics on
OFCstates = decdata.ch_state - decdata.unch_state;

first_OFC_state = nan(ntr,2);

for tr = find(keep_tr)'
    
    % get first OFC state
    idx = find(OFCstates(tr,t0:t_end)~=0,1,'first');
    if ~isempty(idx)
        first_OFC_state(tr,:) = [OFCstates(tr,t0+idx-1), ...
            t_mids(t0+idx-1)];
    end
    
    
    
end



%% show OFC state split by 1st CdN state
figure;
clrs = [0 0 1; 1 0 0];

for s = 1:nsubj
    
    
    subject = subject_names{s};
    
    tr_mis = tr_sampled.(subject).keep_mis;
    ntr_mini = length(tr_mis);
    
    
    counts_mis = hist(first_OFC_state(tr_mis,1),[-1,1]);
    counts_mis = counts_mis/ntr_mini;
    
    counts_match = nan(niters,2);
    for iter = 1:niters
        
        tr_match = tr_sampled.(subject).keep_match(:,iter);
        counts_match(iter,:) = hist(first_OFC_state(tr_match,1),[-1,1]);
        
    end
    
    counts_match = counts_match/ntr_mini;
    
    % plot show proportions
    subplot(2,2,s); hold on
    
    Qs = quantile(counts_match,[.01 .5 .99]);
    
    for i = 1:2
        errorbar(i*2-3, Qs(2,i),Qs(2,i)-Qs(1,i),Qs(3,i)-Qs(2,i),'Color',clrs(i,:));
        scatter(i*2-3,Qs(2,i),'filled','MarkerFaceColor',clrs(i,:))
    end
    scatter([-1.1 .9],counts_mis,'k')
    
    
    ylim([.25 .75])
    
    xlim([-2 2])
    
    set(gca,'XTick',[-1 1],'XTickLabel',{'unchosen','chosen'},'XDir','reverse',...
        'YTick',0:.25:1)
    xlabel('OFC value state')
    ylabel('proportion trials')
    text(1.5,.95,['n=',num2str(ntr_mini)])
    
    title(subject)
    
    % show odds
    subplot(2,2,s+2); hold on
    
    temp = counts_match(:,2)./counts_match(:,1);
    Qs = quantile(temp,[.01 .99]);
    temp_notsig = temp;
    temp_notsig(temp_notsig<=Qs(1) | temp_notsig>=Qs(2)) = NaN;
    
    temp_sig = temp;
    temp_sig(temp_sig>Qs(1) & temp_sig<Qs(2)) = NaN;
    histogram(log2(temp_notsig),'FaceColor',[.7 .7 .7],'BinWidth',.025,'Normalization','probability');
    
    histogram(log2(temp_sig),'FaceColor',[.7 0 .7],'BinWidth',.025,'Normalization','probability');
    y = get(gca,'YLim');
    
    set(gca,'XTick',-2:2,'XTickLabel',{'1 : 4','1 : 2','1 : 1','2 : 1','4 : 1'})
    xlabel('chosen : unchosen')
    
    xlim([-.5 1.5])
    
    b=plot(log2(counts_mis(2)/counts_mis(1))*[1 1],y,'k','LineWidth',2);
    ylabel('norm. trial count')
  
  
    plot([0 0],y,'k:')
  
end


