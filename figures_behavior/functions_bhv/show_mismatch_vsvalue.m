function show_mismatch_vsvalue(bhvdata)
%%
figure;
subject_names = unique(bhvdata.subject);

for s = 1:length(subject_names)
    
   % restrict trials
   tr_keep = strcmp(bhvdata.subject, subject_names{s}) & ... % subject
       bhvdata.trialtype==2 & bhvdata.lever~=0 & ... % free trials
       bhvdata.chose_high==1; % correct only
   
   % get info
   valmax = max(bhvdata.valbin_expval(tr_keep,:),[],2);
   valdiff = abs(diff(bhvdata.valbin_expval(tr_keep,:), [], 2));
   sacc1 = bhvdata.saccloc(tr_keep,1) == bhvdata.lever(tr_keep);
   ntr = length(valmax);
   
   % make vis
   counts = nan(2, 3);
   for val = 1:3%2:4
%        idx = valmax==val;
        idx = valdiff==val;
       counts(:, val) = [sum(sacc1(idx)), sum(~sacc1(idx))]; %counts(:, val-1) 
   end
      
   subplot(1,2,s); 
   bar(counts(2,:)./sum(counts))
   
   xlabel('|\Delta Value|')
   ylabel({'Proportion of not choosing 1st saccade target (correct trials)'})
   
   set(gca,'box','off','YTick',0:.1:.2)

   ylim([0 .2])
  
   % do chi2 test 
   [chi2,df,chi2p] = do_chi2(counts);
   sum(counts)
end
end

