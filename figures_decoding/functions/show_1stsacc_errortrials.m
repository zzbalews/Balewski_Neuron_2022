function show_1stsacc_errortrials(bhvdata, tr_sampled, clrs)
% compare 1st saccade direction between error trials and correct trials, matched for offer
% value
%%
subject_names = fieldnames(tr_sampled);
trtypes = fieldnames(tr_sampled.(subject_names{1}));

correct_side = bhvdata.lever .* bhvdata.chose_high;
% sacc1_match = bhvdata.saccloc(:,1)==bhvdata.lever;
sacc1_match = bhvdata.saccloc(:,1)==correct_side;

figure;
for s = 1:length(subject_names)
    
   subplot(1,2,s); hold on
   
   ll = [];
   for tr = 1:length(trtypes)
       sacc1loc = mean(sacc1_match(tr_sampled.(subject_names{s}).(trtypes{tr})));
       
       ci999 = quantile(sacc1loc, [0.005 0.5 0.995]);
       
       b=bar(tr, ci999(2), 'FaceColor', clrs(tr,:), 'EdgeAlpha', 0, 'FaceAlpha', .7);
       ll(tr) = b;
       e = errorbar(tr, ci999(2), ...
           diff(ci999(1:2)), ...
           diff(ci999(2:3)), 'k.');
       if tr==length(trtypes)
           ll(end+1) = e;
       end
           
   end
   
   plot([0 3],[.5 .5],'k:')
   set(gca,'XTick',[], 'YTick',[0:.5:1])
%    set(gca,'XTick',[1,2],'XTickLabel',{'correct','error'}, 'YTick',[0:.5:1])
   ylim([0 1])
   xlabel('free trials')
   ylabel('Proportion 1^{st} saccade to high value')
   title(subject_names{s})
   
%    legend(ll,{'correct','error','99% CI'},'box','off','Location','northeast')
   
end

end

