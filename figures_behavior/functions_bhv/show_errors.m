function show_errors(bhvdata)
%show how errors are defined

subject_names = unique(bhvdata.subject);
subj_rng = [0.35,0.25];


%% subjval
figure;
for s = 1:length(subject_names)
    
    % restrict to free trials from subject
    subject = subject_names{s};
    tr_subj = strcmp(bhvdata.subject,subject) & ...
        bhvdata.trialtype==2 & ...
        bhvdata.lever~=0 & ~isnan(bhvdata.saccloc(:,1));
    
    
    % get chosen/unchosen value
    lever = bhvdata.lever(tr_subj);
    subjval = bhvdata.subjval_expval(tr_subj,:);
    
    subjval(lever==1,:) = subjval(lever==1,[2,1]); % col1 = ch, col2 = unch
    
    % plot
    valdiff_chosen = subjval(:,1) - subjval(:,2);    
    
    subplot(1,2,s); hold on
    
        ntr_good = valdiff_chosen>=subj_rng(s);
        ntr_meh = valdiff_chosen<subj_rng(s) & ...
            valdiff_chosen>=-subj_rng(s);
        ntr_error = valdiff_chosen<-subj_rng(s);
%     
%     ntr_good = valdiff_chosen>=0;
%     ntr_error = valdiff_chosen<0;
    
        histogram(valdiff_chosen(ntr_good),...
            'BinWidth',.15,'FaceColor',[0 .7 0]);
    
        histogram(valdiff_chosen(ntr_error),...
            'BinWidth',.15,'FaceColor',[.7 0 0]);
    
        histogram(valdiff_chosen(ntr_meh),...
            'BinWidth',.15,'FaceColor',[.3 .3 .3]);
    
    
    y = get(gca,'YLim');
    
    
%     plot([0 0],[y],'k','LineWidth',2);
    
plot(subj_rng(s)*[1 1],[y],'k','LineWidth',2);
plot(-subj_rng(s)*[1 1],[y],'k','LineWidth',2);
    
    
    text(-4,y(2)*.8,[num2str(round(mean(ntr_error)*100,1)),'%'],...
        'Color',[0.7 0 0])
    text(-4,y(2)*.75,['n=',num2str(sum(ntr_error))],...
        'Color',[0.7 0 0])
    
    text(3,y(2)*.8,[num2str(round(mean(ntr_good)*100,1)),'%'],...
        'Color',[0 .7 0])
    text(3,y(2)*.75,['n=',num2str(sum(ntr_good))],...
        'Color',[0 .7 0])
    
    text(-2,y(2)*.8,[num2str(round(mean(ntr_meh)*100,1)),'%'],...
        'Color',[.3 .3 .3])
    text(-2,y(2)*.75,['n=',num2str(sum(ntr_meh))],...
        'Color',[.3 .3 .3])
    
    
    xlim([-5 5])
    
    xlabel('subj.val._{ch} - subj.val._{unch}')
    ylabel('trial count')
    
    title(subject)
    
end

set(gcf,'Position',[50 50 800 300])
end

