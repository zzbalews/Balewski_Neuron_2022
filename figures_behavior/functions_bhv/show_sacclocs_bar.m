function show_sacclocs_bar(bhvdata)
%show 'variable' in 4x4 pic arrangement

%% show subject results (separte plots)
figure;

subject_names = unique(bhvdata.subject);
for s = 1:length(subject_names)
    
    % restrict to ('trialtype') trials from subject
    subject = subject_names{s};
    tr_subj = strcmp(bhvdata.subject,subject) & ...
         bhvdata.trialtype==2 & ...
            bhvdata.lever~=0;
    
        
    sacc_to_lever = [mean(bhvdata.saccloc(tr_subj,1)==bhvdata.lever(tr_subj)),...
        mean(bhvdata.saccloc(tr_subj,2)==bhvdata.lever(tr_subj))];
    sacc_to_opp = [mean(bhvdata.saccloc(tr_subj,1)==-bhvdata.lever(tr_subj)),...
        mean(bhvdata.saccloc(tr_subj,2)==-bhvdata.lever(tr_subj))];
       
    %% plot
    subplot(1,2,s)
    bar([sacc_to_lever; sacc_to_opp]')
    
    title(subject)
    
    set(gca,'XTick',1:2,'XTickLabel',{'1st','2nd'})
    xlabel('saccade')
    ylabel({'proportion trials','saccade==lever direction'})
    
    disp(subject)
    disp(sacc_to_lever)
    
end


end

