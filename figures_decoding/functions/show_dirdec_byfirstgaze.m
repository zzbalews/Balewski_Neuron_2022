function show_dirdec_byfirstgaze(bhvdata,decdata,tr_sampled,decoder,level,tr_sameval)
%%
subject_names = unique(decdata.subject);
nsubj = length(subject_names);


% useful vars
t_mids = decdata.t_mids;
nmids = length(t_mids);
XLIM = [-500 1000];
niters = size(tr_sampled.(subject_names{1}).keep_match,2);

if strcmp(decoder,'direction')
    
    
    % relabel by chosen direction
    decdata.postprob_ch = decdata.postprob(:,:,1);
    idxR = bhvdata.lever==1;
    decdata.postprob_ch(idxR,:) = decdata.postprob(idxR,:,2);
    
    
    pp{1} = decdata.postprob_ch;
    YLIM = [.3 .8];%[.4 .75];
    YTICK = 0:.1:1;
    chance = .5;
    ROWS = 1;
    
    yy = {'chosen'};
    metric = 'posterior probability';
    
    clrs = [1 .6 0;.2 .2 .2];
    
elseif strcmp(decoder,'value')
    
    if strcmp(level,'chosen')
    pp{1} = decdata.ch_ppd;
    yy = {'chosen'};
    clrs = [1 0 0;.2 .2 .2];
    elseif strcmp(level,'unchosen')
    pp{1} = decdata.unch_ppd;
    yy = {'unchosen'};
    clrs = [0 0 1; .2 .2 .2];
    end
    YLIM = [-75 200];
    YTICK = 0:100:200;
    chance = 0;
    ROWS = 1;
    
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
        
        % match trials
        temp = nan(niters,nmids);
        for j = 1:nmids
            HH = reshape(tr_sampled.(subject).keep_match,1,[]);
            KK = pp{r}(HH,j)';
            LL = reshape(KK,[],niters);
            temp(:,j) = mean(LL);
        end

        Y = quantile(temp,[.0005,.5,.9995],1);
        shadedErrorBar(t_mids, Y(2,:),...
            [Y(1,:)-Y(2,:);Y(2,:)-Y(3,:)],'lineprops',{'Color',clrs(1,:)});
      
        % add match RT:
        rt = median(bhvdata.rt(tr_sampled.(subject).keep_match));
        rt999 = quantile(rt, [0.0005 0.5 0.9995]);
%         e = errorbar(rt999(2), .7, [], [], ...
%            diff(rt999(1:2)), ...
%            diff(rt999(2:3)), ...
%            'o','MarkerSize',2,'MarkerFaceColor',clrs(1,:),'Color',clrs(1,:));
       
        % mis trials
        temp = nan(niters,nmids);
        for j = 1:nmids
            HH = reshape(tr_sampled.(subject).keep_mis,1,[]);
            KK = pp{r}(HH,j)';
            LL = reshape(KK,[],niters);
            temp(:,j) = mean(LL);
        end

        Y = quantile(temp,[.0005,.5,.9995],1);
        shadedErrorBar(t_mids, Y(2,:),...
            [Y(1,:)-Y(2,:);Y(2,:)-Y(3,:)],'lineprops',{'Color',clrs(2,:)});
        
        % add mis RT:
        rt = median(bhvdata.rt(tr_sampled.(subject).keep_mis));
        rt999 = quantile(rt, [0.0005 0.5 0.9995]);
%         e = errorbar(rt999(2), .7, [], [], ...
%            diff(rt999(1:2)), ...
%            diff(rt999(2:3)), ...
%            'o','MarkerSize',2,'MarkerFaceColor',clrs(2,:),'Color',clrs(2,:)); 
       
        
        a = plot([NaN NaN],[.5 .5],'Color',clrs(1,:));
        b = plot([NaN NaN],[.5 .5],'Color',clrs(2,:));

%         b = plot(t_mids, mean(pp{r}(tr_sampled.(subject).keep_mis,:)),'Color',clrs(2,:),'LineWidth',1);
        
        % add sacc refs
        scatter(refs,(YLIM(2)-scale)*[1 1 1],'vk','filled')
        text(refs,YLIM(2)*[1 1 1],{'1st','2nd','RT'},'HorizontalAlignment','center')
% % % 
% % % refs_y = [YLIM(1) YLIM(2)-2*scale]';
% % %         plot(refs.*[1;1],repmat(refs_y,1,3),'k--')

        % format
        
        text(-250,YLIM(2)*0.9,['n=',num2str(size(tr_sampled.(subject).keep_mis,1))],'HorizontalAlignment','center')
%         ll = legend([a,b],{'lever','other'},'box','off');
%         title(ll,'1st sacc. direction')
        xlim(XLIM)
        ylim(YLIM)
        
        set(gca,'XTick',-500:500:1500,'XTickLabel',-.5:.5:1.5,...
            'YTick',YTICK)

        xlabel('Time from pictures (s)')
        ylabel([yy{r},' ',decoder,' decoder: ',metric])
    end
end


end

