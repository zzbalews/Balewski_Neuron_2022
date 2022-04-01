function show_1stsacc_errortrials(bhvdata, dirdata, tr_sampled, clrs)
% compare CN dir dec post prob between error trials and correct trials, matched for offer
% value
%%
subject_names = fieldnames(tr_sampled);
trtypes = fieldnames(tr_sampled.(subject_names{1}));

nmids = size(dirdata.postprob_ch,2);
niters = size(tr_sampled.(subject_names{1}).(trtypes{1}),2);

t_mids = dirdata.t_mids;

figure;
for s = 1:length(subject_names)
    
   subplot(1,2,s); hold on
       plot(t_mids([1 end]),[.5 .5],'k:')

    plot([0 0],[0 1],'k:')
   ll = [];
   
   for tr = 1:length(trtypes)
       tic
       postprob = nan(niters,nmids);
       for m = 1:nmids % work through each time point
           idx = reshape(tr_sampled.(subject_names{s}).(trtypes{tr}), 1, []);
           recover = reshape(dirdata.postprob_ch(idx,m), [], niters);
           postprob(:,m) = mean(recover);
       end
          
      
       
       ci999 = quantile(postprob, [0.005 0.5 0.995],1);
       
       if strcmp(trtypes{tr}, 'keep_bad')
           ci999 = 1 + -1 * ci999;
       end
      
       b = shadedErrorBar(t_mids, ci999(2,:), ...
           [diff(ci999(1:2,:));diff(ci999(2:3,:))],...
           'lineprops',{'Color',clrs(tr,:)});
       
       ll(tr) = b.mainLine;
       
  toc
   end
   
   ylim([.35 .75])
   xlim([-500 1000])
   xlabel('free trials')
   ylabel('Direction decoder: high value post prob')
   title(subject_names{s})
   
   legend(ll,{'correct','error'},'box','off','Location','northeast')
   
end

end

