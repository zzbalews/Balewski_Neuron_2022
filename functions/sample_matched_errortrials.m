function output = sample_matched_errortrials(bhvdata, tr_error, tr_correct, niters, ifshow)
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
    tr_subj = strcmp(bhvdata.subject,subject);
    
    % split by error/correct
    tr_bad = find(tr_subj & tr_error); % trial#
    tr_good = find(tr_subj & tr_correct);% trial#
    
    disp(['# bad=',num2str(length(tr_bad))])
    disp(['# good=',num2str(length(tr_good))])
    
    % possible picture configurations
    freetypes = [repelem([1;2;3;4],4),repmat([1;2;3;4],4,1)];
    
    keep_bad = nan(0,niters);
    keep_good = nan(0,niters);
    
    
    for u = 1:size(freetypes,1)
        
        % pull idx, within 'tr_good' or 'tr_bad'
        idx_bad = sum(bhvdata.valbin_expval(tr_bad,:)==freetypes(u,:),2)==2;
        idx_good = sum(bhvdata.valbin_expval(tr_good,:)==freetypes(u,:),2)==2;
        
        % use same # match trials as mis
        N = min(sum(idx_bad),sum(idx_good));
        
        % update bad trials: sample with replacement
        temp = tr_bad(find(idx_bad));
        blk = nan(N,niters);
        for i = 1:niters
            add_these=temp(randsample(sum(idx_bad),N,'true'));
            blk(:,i) = add_these;
        end
        keep_bad = cat(1,keep_bad, blk);
        
        
        % update good trials: sample with replacement
        temp = tr_good(find(idx_good));
        blk = nan(N,niters);
        for i = 1:niters
            add_these=temp(randsample(sum(idx_good),N,'true'));
            blk(:,i) =add_these;
        end
        keep_good = cat(1,keep_good,blk);
        
    end
    
    % update output
    output.(subject).keep_good = keep_good;
    output.(subject).keep_bad = keep_bad;
    
%     if ifshow
%         %         subplot(1,2,s); hold on
%         %
%         %         histogram(bhvdata.rt(keep_match,:),'FaceColor',clrs(1,:),'Normalization','probability','BinWidth',50);
%         %         histogram(bhvdata.rt(keep_mis,:),'FaceColor',clrs(2,:),'Normalization','probability','BinWidth',50);
%         %
%         %         med_match = mean(mean(bhvdata.rt(keep_match)));
%         %         med_mis = mean(bhvdata.rt(keep_mis));
%         %
%         %         title([subject,': \Delta=',num2str(round(med_mis - med_match)),'ms'])
%         %
%         %         ll = legend({'lever','other'},'box','off');
%         %         title(ll,'1st sacc. direction')
%         %         xlim([0 1500])
%         %
%         %         xlabel('RT (ms)')
%         %         ylabel('norm. trial count')
%         
%         
%         subplot(1,2,s); hold on
%         m_good = nan(niters,1);
%         for i = 1:niters
%             m_good(i) = mean(bhvdata.rt(keep_good(:,i)));
%         end
%         histogram(m_good,'EdgeAlpha',0,'FaceColor',[.2 .2 .2],'Normalization','probability')
%         
%         Y = get(gca,'YLim');
%         
%         m_bad = mean(bhvdata.rt(keep_bad));
% %         plot(m_bad*[1 1],Y,'Color',[1 0 0],'LineWidth',2)
%         histogram(m_bad,'EdgeAlpha',0,'FaceColor',[1 0 0],'Normalization','probability')
% 
%         
%         X = get(gca,'XLim');
%         
%         set(gca,'XTick',500:100:1500,'XTickLabel',.5:.1:1.5,...
%             'YTick',linspace(Y(1),Y(2),3))
%         xlabel('response time (s)')
%         ylabel('norm. trial count')
%         ylim(Y)
%         xlim(X)
%         text(mean(X),mean(Y),['\Delta=',num2str(round(mean(m_bad)-mean(m_good))),'ms'],'HorizontalAlignment','left')
%     end
end

end

