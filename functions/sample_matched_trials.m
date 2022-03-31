function output = sample_matched_trials(bhvdata,keep_tr,niters,ifshow)
%%

if ifshow
    figure;
end

subject_names = unique(bhvdata.subject);
clrs = [0 .6 0; .7 0 .7];

output = struct();

for s = 1:length(subject_names)
    
    % restrict trials
    subject = subject_names{s};
    tr_subj = strcmp(bhvdata.subject,subject) & ...
        bhvdata.trialtype==2 & ...
        bhvdata.lever~=0 & ~isnan(bhvdata.saccloc(:,1)) & ...
        keep_tr;
    
    % split by 1st saccade
    tr_match = find(tr_subj & bhvdata.lever==bhvdata.saccloc(:,1)); % trial#
    tr_mis = find(tr_subj & -bhvdata.lever==bhvdata.saccloc(:,1));% trial#
    
% %     % split by first CdN state = chosen vs not
% %     tr_match = find(tr_subj & bhvdata.first_CdN_chosen==1); % trial#
% %     tr_mis = find(tr_subj & ~bhvdata.first_CdN_chosen==1);% trial#
% %     
% %     % split by fist OFC state = chosen vs unchosen
% %     tr_match = find(tr_subj & bhvdata.first_OFC==1); % trial#
% %     tr_mis = find(tr_subj & bhvdata.first_OFC==-1);% trial#

    disp(['#match=',num2str(length(tr_match))])
    disp(['#mis=',num2str(length(tr_mis))])
    
    % possible picture configurations
    freetypes = [repelem([1;2;3;4],4),repmat([1;2;3;4],4,1)];
    
    % use all mismatch trials
    %     keep_mis = tr_mis;
    keep_mis = nan(0,niters);
    keep_match = nan(0,niters);
    
    
    for u = 1:size(freetypes,1)
        
        % pull idx, within 'tr_match' or 'tr_mis'
        idx_match = sum(bhvdata.valbin_expval(tr_match,:)==freetypes(u,:),2)==2;
        idx_mis = sum(bhvdata.valbin_expval(tr_mis,:)==freetypes(u,:),2)==2;
        
        
        % use same # match trials as mis
        N = min(sum(idx_mis),sum(idx_match));
        
        % update mis trials
        temp = tr_mis(find(idx_mis));
%         add_these=temp(randsample(sum(idx_mis),N,'true'));
%         keep_mis = cat(1,keep_mis,add_these);
blk = nan(N,niters);
for i = 1:niters
    add_these=temp(randsample(sum(idx_mis),N,'true'));
    blk(:,i) = add_these;
end
keep_mis = cat(1,keep_mis, blk);
        
        
        % update match trials
        temp = tr_match(find(idx_match));
        blk = nan(N,niters);
        for i = 1:niters
            add_these=temp(randsample(sum(idx_match),N,'true'));
            blk(:,i) =add_these;
        end
        keep_match = cat(1,keep_match,blk);
        
    end
    
    % update output
    output.(subject).keep_match = keep_match;
    output.(subject).keep_mis = keep_mis;
    
    if ifshow
%         subplot(1,2,s); hold on
%         
%         histogram(bhvdata.rt(keep_match,:),'FaceColor',clrs(1,:),'Normalization','probability','BinWidth',50);
%         histogram(bhvdata.rt(keep_mis,:),'FaceColor',clrs(2,:),'Normalization','probability','BinWidth',50);
%         
%         med_match = mean(mean(bhvdata.rt(keep_match)));
%         med_mis = mean(bhvdata.rt(keep_mis));
%         
%         title([subject,': \Delta=',num2str(round(med_mis - med_match)),'ms'])
%         
%         ll = legend({'lever','other'},'box','off');
%         title(ll,'1st sacc. direction')
%         xlim([0 1500])
%         
%         xlabel('RT (ms)')
%         ylabel('norm. trial count')
        
        
        subplot(1,2,s); hold on
        m_match = nan(niters,1);
        for i = 1:niters
            m_match(i) = mean(bhvdata.rt(keep_match(:,i)));
        end
        histogram(m_match,'EdgeAlpha',0,'FaceColor',[.2 .2 .2],'Normalization','probability')
        
        Y = get(gca,'YLim');
        
        m_mis = mean(bhvdata.rt(keep_mis));
        plot(m_mis*[1 1],Y,'Color',[1 0 0],'LineWidth',2)
        
                X = get(gca,'XLim');

        set(gca,'XTick',500:100:1500,'XTickLabel',.5:.1:1.5,...
            'YTick',linspace(Y(1),Y(2),3))
        xlabel('response time (s)')
        ylabel('norm. trial count')
        ylim(Y)
        xlim(X)
        text(mean(X),mean(Y),['\Delta=',num2str(round(m_mis-mean(m_match))),'ms'],'HorizontalAlignment','left')
    end
end

end

