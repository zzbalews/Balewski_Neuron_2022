function show_responsetime_hists(bhvdata)
%show 'variable' histogram across all trials

%% set colors/style
clr_free = [0 .4 .8];
clr_forced = [0 0 0];
clr_sacc1 = [.6 0 .6];
clr_sacc2 = [1 .6 0];

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
    
    track_sacc1time = bhvdata.MLsacc(tr_subj_free,1);
    track_sacc2time = bhvdata.MLsacc(tr_subj_free,2);
    
    subplot(1,2,s)
    hold on;
    rt_bin = 50;
    sacc_bin = 10;
    
    histogram(track_RT_forced,'Binwidth',rt_bin,...
        'Normalization','probability',...
        'FaceColor',clr_forced,'EdgeAlpha',0,'FaceAlpha',.5);%'DisplayStyle','stairs','EdgeColor',clr_forced,'EdgeAlpha',.7);
    histogram(track_RT_free,'Binwidth',rt_bin,...
        'Normalization','probability',...
        'FaceColor',clr_free,'EdgeAlpha',0,'FaceAlpha',.5);% 'DisplayStyle','stairs','EdgeColor',clr_free,'EdgeAlpha',.7);
    
    histogram(track_sacc1time,'BinWidth',sacc_bin,...
        'Normalization','probability',...
        'FaceColor',clr_sacc1,'EdgeAlpha',0,'FaceAlpha',.9);
    histogram(track_sacc2time,'BinWidth',sacc_bin,...
        'Normalization','probability',...
        'FaceColor',clr_sacc2,'EdgeAlpha',0,'FaceAlpha',.9);
    
%     title(subject)
    
    xlim([0 1500])
    ylim([0 .22])
    
    xlabel('time from pics on (ms)')
    ylabel('normalized trial count')
    
    legend({'rt forced','rt free',...
        'sacc1 free',...
        'sacc2 free'},...
        'box','off')
end



% %% show prop trials where saccade = chosen dir
% figure;
% 
% for s = 1:length(subject_names)
%     
%     % restrict to ('trialtype') trials from subject
%     subject = subject_names{s};
%     tr_subj_free = strcmp(bhvdata.subject,subject) & ...
%         bhvdata.trialtype==2 & ...
%         bhvdata.lever~=0;
%     
%     ntr = sum(tr_subj_free);
%     
%     % prop trials where saccades match lever?
%     
%     Nsacc1_ch = sum(tr_subj_free & bhvdata.saccloc_chosen(:,1)==1);
%     Nsacc1_unch = sum(tr_subj_free & bhvdata.saccloc_chosen(:,1)==0);
%     
%     Nsacc2_ch = sum(tr_subj_free & bhvdata.saccloc_chosen(:,2)==1);
%     Nsacc2_unch = sum(tr_subj_free & bhvdata.saccloc_chosen(:,2)==0);
%     
%     
%     vars = [Nsacc1_ch,Nsacc1_unch,Nsacc2_ch,Nsacc2_unch]/ntr;
%     
%     x = [.8 1.2 1.8 2.2];
%     scale = [1 1 1 1];
%     clrs = min(max([clr_sacc1; clr_sacc1; clr_sacc2; clr_sacc2] - 0.2*repmat([-1;1;-1;1],1,3),0),1);
%     
%     subplot(1,2,s); hold on
%     
%     for b = 1:4
%         bar(x(b),scale(b)*vars(b),.4,'FaceColor',clrs(b,:),'EdgeAlpha',0)
%     end
%     ylim([0 1])
%     
%     set(gca,'XTick',1:2,'XTickLabel',{'1^{st}','2^{nd}'},'box','off');
%     
%     xlabel('saccade after pics on')
%     ylabel('proportion trials')
%     
%     legend({'to chosen','opposite','to chosen','opposite'},'box','off')
%     
%     sum(vars(1:2))
%     sum(vars(3:4))
% end


end
