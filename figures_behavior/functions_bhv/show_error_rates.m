function show_error_rates(bhvdata)
% show error rate, usng subjective value and value bins

%% show subject results (separate plots)
% figure;

subject_names = unique(bhvdata.subject);

  % check if picked higher subjective value
   temp = bhvdata.subjval_expval;
   temp(bhvdata.lever==1,:) = temp(bhvdata.lever==1,[2,1]);
   
   bhvdata.chose_high_subjval = temp(:,1)>temp(:,2);
   
  
   track_err = nan(2,2); % subjects x [subjval binval]
for s = 1:length(subject_names)
    
   % restrict to free trials
   subject = subject_names{s};
   
   tr_subj_free = strcmp(bhvdata.subject,subject) & ...
       bhvdata.trialtype==2 & ...
       bhvdata.lever~=0;
   
   % error rate using subjective value
   track_err(s,1) = mean(bhvdata.chose_high_subjval(tr_subj_free)==0);
   % error rate using value bin
    track_err(s,2) = mean(bhvdata.chose_high(tr_subj_free)==-1);          
    
end

round(track_err*100,1)
end

