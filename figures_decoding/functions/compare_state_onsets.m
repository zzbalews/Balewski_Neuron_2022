function compare_state_onsets(bhvdata,dirdata_pics,dirdata_choice,keep_tr)

subject_names = unique(dirdata_pics.subject);
nsubj = length(subject_names);

t_mids = dirdata_pics.t_mids;
nmids = length(t_mids);

t_pics = find(t_mids>=0 & t_mids<400);
t_choice = find(t_mids>=-400 & t_mids<0);

figure;
for S = 1:nsubj
    
    subject = subject_names{S};
    idx_subj = find(strcmp(bhvdata.subject,subject));
    
   ntr = length(idx_subj);
   
     pics_states = dirdata_pics.postprob_ch_disc(idx_subj,:);
   choice_states = dirdata_choice.postprob_ch_disc(idx_subj,:);
   
   pics_state_onset = nan(ntr,1);
   choice_state_onset = nan(ntr,1);
   
   for tr = 1:ntr
      
       temp = pics_states(tr,t_pics);
       onset = find(abs(temp),1);
       if ~isempty(onset)
       pics_state_onset(tr) = t_mids(onset + t_pics(1) -1);
       end
       
       temp = choice_states(tr,t_choice);
       onset = find(abs(temp),1);
       if ~isempty(onset)
       choice_state_onset(tr) = t_mids(onset + t_pics(1) -1);
       end
   end

   subplot(1,2,S); hold on
   a=histogram(choice_state_onset,...
       'Normalization','probability','BinWidth',20,...
       'FaceColor',[.1 .1 .1]);
   b=histogram(pics_state_onset,...
       'Normalization','probability','BinWidth',20,...
       'FaceColor',[0 .4 .8]);
   
   Y = get(gca,'YLim');
   set(gca,'XTick',0:200:400,'XTickLabel',0:.2:.4,'YTick',linspace(0,round(Y(2),1),3))
   xlabel('first state onset (s)')
   ylabel('norm. trial count')
   ll = legend([b,a],{'pics','choice'},'box','off');
   title(ll,'sync trials to:')
end


end

