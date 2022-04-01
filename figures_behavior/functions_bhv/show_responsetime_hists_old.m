function show_responsetime_hists(bhvdata)
%show 'variable' histogram across all trials

%% set colors/style
clr_free = [0 .4 .8];
clr_forced = [0 0 0];

%% do saccs match chosen dir?
bhvdata.saccloc_chosen = double(bhvdata.saccloc==bhvdata.lever);
bhvdata.saccloc_chosen(isnan(bhvdata.saccloc))=NaN;
%% show subject results (separte plots)
figure;

subject_names = unique(bhvdata.subject);
for s = 1:length(subject_names)
    
    % restrict to ('trialtype') trials from subject
    subject = subject_names{s};
    tr_subj_free = strcmp(bhvdata.subject,subject) & ...
        bhvdata.trialtype==2 & ...
        bhvdata.lever~=0;
    
    tr_subj_forced = strcmp(bhvdata.subject,subject) & ...
        bhvdata.trialtype==1 & ...
        bhvdata.lever~=0;
    
    
    
    track_RT_free = bhvdata.rt(tr_subj_free);
    track_RT_forced = bhvdata.rt(tr_subj_forced);
    
    
    track_sacc1time_free = bhvdata.MLsacc(tr_subj_free & bhvdata.saccloc_chosen(:,1)==1,1);
    track_sacc1time_forced = bhvdata.MLsacc(tr_subj_forced & bhvdata.saccloc_chosen(:,1)==1,1);
    
    track_sacc2time_free = bhvdata.MLsacc(tr_subj_free,2);
    track_sacc2time_forced = bhvdata.MLsacc(tr_subj_forced,2);
    
    
    subplot(1,2,s)
    hold on;
    rt_bin = 60;
     histogram(track_RT_forced,'Binwidth',rt_bin,...
        'Normalization','probability','FaceColor',clr_forced,...
        'EdgeAlpha',0);
    histogram(track_RT_free,'Binwidth',rt_bin,...
        'Normalization','probability','FaceColor',clr_free,...
        'EdgeAlpha',0);
   
    
    sacc_bin = 15;
    histogram(track_sacc1time_forced,'Binwidth',sacc_bin,...
        'Normalization','probability','EdgeColor',clr_forced,...
        'EdgeAlpha',.7,'DisplayStyle','stairs','LineWidth',1);
    
    histogram(track_sacc1time_free,'Binwidth',sacc_bin,...
        'Normalization','probability','EdgeColor',clr_free,...
        'EdgeAlpha',.7,'DisplayStyle','stairs','LineWidth',1);
%     
%     histogram(track_sacc2time_forced,'Binwidth',sacc_bin,...
%         'Normalization','probability','FaceColor',clr_forced,...
%         'EdgeAlpha',0,'FaceAlpha',.2);
%     
%     histogram(track_sacc2time_free,'Binwidth',sacc_bin,...
%         'Normalization','probability','FaceColor',clr_free,...
%         'EdgeAlpha',0,'FaceAlpha',.2);
    
    title(subject)
    
    xlim([0 1750])
    ylim([0 .28])
    
    xlabel('time from pics on (ms)')
    ylabel('normalized trial count')
    
    legend({'rt forced','rt free',...
        'sacc1 forced','sacc1 free'},...'sacc2 forced','sacc2 free'},...
        'box','off')
end


%%
%% show subject results (separte plots)
figure;

subject_names = unique(bhvdata.subject);
for s = 1:length(subject_names)
    
    % restrict to ('trialtype') trials from subject
    subject = subject_names{s};
    tr_subj_free = strcmp(bhvdata.subject,subject) & ...
        bhvdata.trialtype==2 & ...
        bhvdata.lever~=0;
    
    tr_subj_forced = strcmp(bhvdata.subject,subject) & ...
        bhvdata.trialtype==1 & ...
        bhvdata.lever~=0;
    
   
    track_sacc1time_free = bhvdata.MLsacc(tr_subj_free & bhvdata.saccloc_chosen(:,1)==1,1);
    track_sacc1time_forced = bhvdata.MLsacc(tr_subj_forced & bhvdata.saccloc_chosen(:,1)==1,1);
    
    track_sacc1time_free_unch = bhvdata.MLsacc(tr_subj_free & bhvdata.saccloc_chosen(:,1)==0,1);
    
    track_sacc2time_free = bhvdata.MLsacc(tr_subj_free & bhvdata.saccloc_chosen(:,2)==1,2);
    track_sacc2time_forced = bhvdata.MLsacc(tr_subj_forced & bhvdata.saccloc_chosen(:,2)==1,2);
    
    track_sacc2time_free_unch = bhvdata.MLsacc(tr_subj_free & bhvdata.saccloc_chosen(:,2)==0,2);
    
    track_RT_free = bhvdata.rt(tr_subj_free);
    track_RT_forced = bhvdata.rt(tr_subj_forced);
    
    vars = {track_RT_free,track_sacc1time_free,track_sacc1time_free_unch,...
       track_sacc2time_free, track_sacc2time_free_unch, };
    scale = [1 1 -1 1 -1];
    width = 30*[1 1 1 1 1];
    alpha = [1 .7 .7 .9 .9];
    
    subplot(1,2,s);hold on
    for v = 1:length(vars)
    [N,X] = hist(vars{v},0:width(v):2000);
    bar(X,scale(v)*N,1,'FaceAlpha',alpha(v));
    end
    xlim([0 1500])
    
    legend({'sacc1','sacc1 opp','sacc2','sacc2 opp','rt'},'box','off')
end

end
