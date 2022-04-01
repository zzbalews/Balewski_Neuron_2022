function show_dirdec_byfirstgaze_mini(bhvdata,decdata,tr_sampled,decoder,tr_sameval)
%%
subject_names = unique(decdata.subject);
nsubj = length(subject_names);


% useful vars
t_mids = decdata.t_mids;
nmids = length(t_mids);
XLIM = [-500 1000];
niters = size(tr_sampled.(subject_names{1}).keep_match,2);
% clrs = [1 .6 0;.2 .2 .2];
% clrs = [1 0 0;.2 .2 .2];
clrs = [1 0 0 ; 0 0 1];

if strcmp(decoder,'direction')
    
    
    % relabel by chosen direction
    decdata.postprob_ch = decdata.postprob(:,:,1);
    idxR = bhvdata.lever==1;
    decdata.postprob_ch(idxR,:) = decdata.postprob(idxR,:,2);
    
    
    pp{1} = decdata.postprob_ch;
    YLIM = [.45 .85];
    YTICK = 0:.1:1;
    chance = .5;
    ROWS = 1;
    
    yy = {'chosen'};
    metric = 'posterior probability';
elseif strcmp(decoder,'value')
    pp{1} = decdata.ch_ppd;
    pp{2} = decdata.unch_ppd;
    YLIM = [-50 200];
    YTICK = 0:100:200;
    chance = 0;
    ROWS = 1;
    
    yy = {'chosen','unchosen'};
    metric = '%\Delta posterior probability';
end

scale = diff(YLIM)*.05;

%% get median sacc & RT times
    
figure;
for S = 1:nsubj
    subject = subject_names{S};
    idx = tr_sameval & strcmp(bhvdata.subject,subject);
    refs = [nanmedian(bhvdata.MLsacc(idx,1:2)),...
        nanmedian(bhvdata.rt(idx))];
    
    for r = 1:ROWS
        % plot avg post prob
        subplot(ROWS,2,S+2*(r-1)); hold on
        plot(XLIM,chance*[1 1],'k:')
        plot([0 0],YLIM,'k:')
        
        match = pp{r}(tr_sampled.(subject).keep_match,:);
        mis = pp{r}(tr_sampled.(subject).keep_mis,:);
        ntr_mini = size(match,1);
        
        a=shadedErrorBar(t_mids, mean(match), std(match)/sqrt(ntr_mini),...
            'lineprops',{'Color',clrs(1,:)});
        b=shadedErrorBar(t_mids, mean(mis), std(mis)/sqrt(ntr_mini),...
            'lineprops',{'Color',clrs(2,:)});
        % add sacc refs
        scatter(refs,(YLIM(2)-scale)*[1 1 1],'vk','filled')
        text(refs,YLIM(2)*[1 1 1],{'1st','2nd','RT'},'HorizontalAlignment','center')

        refs_y = [YLIM(1) YLIM(2)-2*scale]';
%         plot(refs.*[1;1],repmat(refs_y,1,3),'k--')

        % format
        
        text(-250,mean(YLIM),['n=',num2str(length(tr_sampled.(subject).keep_mis))],'HorizontalAlignment','center')
        if S==1
        ll = legend([a.mainLine,b.mainLine],{'chosen','unchosen'},'box','off','Location','northwest');
        title(ll,'1st OFC value state')
        end
        xlim(XLIM)
        ylim(YLIM)
        
        set(gca,'XTick',-500:500:1500,'XTickLabel',-.5:.5:1.5,...
            'YTick',YTICK)

        xlabel('Time from pictures (s)')
        ylabel([yy{r},' ',decoder,' decoder: ',metric])
    end
end


end

