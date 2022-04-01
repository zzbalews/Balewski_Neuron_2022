function compare_LOO_peaks(LOO_pics,LOO_choice)

avgacc_pics = LOO_pics.acc;
avgacc_choice = LOO_choice.acc;
% avgacc_pics = LOO_pics.acc_quants(:,:,4);
% avgacc_choice = LOO_choice.acc_quants(:,:,4);

 [~,peak_pics] = max(avgacc_pics, [],2);
 [~,peak_choice] = max(LOO_choice.acc, [],2);
 
 nses = size(avgacc_pics,1);
 acc_pics = nan(nses,1);
 acc_choice = nan(nses,1);
 
 for i = 1:nses
    
     acc_pics(i) = avgacc_pics(i,peak_pics(i));
     acc_choice(i) = avgacc_choice(i,peak_choice(i));   
     
 end
 
subject_names = unique(LOO_pics.subject);
nsubj = length(subject_names);

for S = 1:nsubj
    
    subject = subject_names{S};
    idx = strcmp(LOO_pics.subject,subject);

    mean(acc_pics(idx));
    mean(acc_choice(idx));
    
    [h,p,ci,stats] = ttest(acc_pics(idx),acc_choice(idx));
    disp([subject,': accuracy @ pics vs @ choice'])
    disp(['   t(',num2str(stats.df),')=',num2str(stats.tstat),',p=',num2str(p)])

end

%%



avgaccQ_pics = LOO_pics.acc_quants;
avgaccQ_choice = LOO_choice.acc_quants;


 accQ_pics = nan(nses,4);
 accQ_choice = nan(nses,4);
 t_mids = LOO_pics.t_mids;
 
 clrs = [0 0 0; 1 0 0 ; .5 0 1; 0 0 1];
 figure; hold on;
      plot([0 1],[0 1],'k:','LineWidth',1)

 for i = 1:nses
     t_p = peak_pics(i);    
    t_c = peak_choice(i);
    
    if strcmp(LOO_pics.subject(i),'Chap')
        clr = [1 0 1];
    else
        clr = [0 .7 0];
    end
        scatter(avgacc_pics(i,t_p), avgacc_choice(i,t_c),'filled','MarkerFaceAlpha',.4,'MarkerFaceColor',clr)
%     scatter(avgaccQ_pics(i,t_p,4), avgaccQ_choice(i,t_c,4),'filled','MarkerFaceAlpha',.4,'MarkerFaceColor',clr)

 end
 
 a=scatter([NaN],[1],'filled','MarkerFaceAlpha',.5,'MarkerFaceColor',[1 0 1]);
 b=scatter([NaN],[1],'filled','MarkerFaceAlpha',.5,'MarkerFaceColor',[0 .7 0]);
 
     legend([a,b],{'Chap','George'},'box','off','Location','northwest')
     daspect([1 1 1])
     xlim([.4 1])
     ylim([.4 1])
     
     xlabel('peak LOO accuracy @ pics')
     ylabel('peak LOO accuracy @ choice')
     
     set(gca,'XTick',.5:.25:1,'YTick',.5:.25:1)
set(gcf,'Position',[100 100 250 250])

% % % %%
% % % track_LOOacc = nan(nses,1);
% % % for i = 1:nses
% % %  t_p = peak_pics(i);    
% % %     track_LOOacc(i) = avgacc_pics(i,t_p);
% % % end
% % % figure; hold on
% % % idx = strcmp(LOO_pics.subject,'George');
% % % scatter(ones(sum(idx),1),track_LOOacc(idx),'filled',...
% % % 'jitter','on',...
% % %     'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',.4)
% % % 
% % % plot([0 2],[.5 .5],'k:')
% % % ylim([.4 1])
% % % xlim([.5 1.5])
% % % ylabel({'Direction decoder','leave-one-out accuracy'})
% % % set(gca,'XTick',1,'XTickLabel','subject G',...
% % %     'YTick',0:.25:1)
% % % 
% % % set(gcf,'Position',[100 100 150 300])
end
