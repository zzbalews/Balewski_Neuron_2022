function show_avg_postprob(bhvdata,dirdata,keep_tr)
%%
subject_names = unique(dirdata.subject);
nsubj = length(subject_names);

% relabel by chosen direction
dirdata.postprob_ch = dirdata.postprob(:,:,1);
idxR = bhvdata.lever==1;
dirdata.postprob_ch(idxR,:) = dirdata.postprob(idxR,:,2);

at_ch = bhvdata.saccloc==bhvdata.lever;
at_unch = bhvdata.saccloc==-bhvdata.lever;
bhvdata.saccloc_target = nan(size(bhvdata.saccloc));
bhvdata.saccloc_target(at_ch) = 1;
bhvdata.saccloc_target(at_unch) = -1;

% useful vars
t_mids = dirdata.t_mids;
XLIM = [-500 1200];

figure;
for S = 1:nsubj
    % restrict trials
    subject = subject_names{S};
    idx = find(strcmp(bhvdata.subject,subject_names{S}) & ...
        keep_tr);
    
    ntr = length(idx);
    
    % plot avg post prob
    subplot(2,2,S); hold on
    plot(t_mids, mean(dirdata.postprob_ch(idx,:)),'LineWidth',2,'Color',[.7 0 .7]);
    
    plot([-2000 2000],[.5 .5],'k:')
    
    xlim(XLIM)
    ylim([.45 .75])
    
    % get current gaze location
    gaze_ch = sum(bhvdata.location_target(idx,:)==bhvdata.lever(idx))';
    gaze_unch = sum(bhvdata.location_target(idx,:)==-bhvdata.lever(idx))';
    gaze_na = sum(bhvdata.location_target(idx,:)==0)';
    
    % plot gaze location
    subplot(2,2,S+2); hold on
    %         A = area((-2000:2000)',[gaze_ch,gaze_unch]/ntr);
    %         A(1).FaceColor = [0 .2 .7];
    %         A(2).FaceColor = [.7 .2 0];
    plot(-2000:2000,gaze_ch/ntr,'LineWidth',2,'Color',[0 .2 .7]);
    plot(-2000:2000,gaze_unch/ntr,'LineWidth',2,'Color',[.7 .2 0]);
    xlim(XLIM)
    ylim([0,1])
    
    % add plot refs
    rt_y = [.7,.9];
    for i = 1:2
        subplot(2,2,S+2*(i-1)); hold on
        
        % pics
        plot([0 0],[0 1],'k:')
        
        % saccs
        sacc1 = nanmedian(bhvdata.MLsacc(idx,1));
        sacc2 = nanmedian(bhvdata.MLsacc(idx,2));
        
        plot(sacc1*[1 1],[0 1],'Color',[.5 .5 .5])
        plot(sacc2*[1 1],[0 1],'Color',[.5 .5 .5])
        
        % rt
        medrt = nanmedian(bhvdata.rt(idx));
        scatter(medrt,rt_y(i),'kv','filled')
        if i==1
            text(medrt,rt_y(i)+.025,num2str(round(medrt)),...
                'HorizontalAlignment','center')
        end
        
        % labels
        xlabel('time from pics (ms)')
        if i==1
            ylabel({'direction decoder','posterior probability'})
        elseif i==2
            ylabel({'gaze location','proportion trials'})
        end
        title(subject)
    end
    
end

%% find avg postprob around saccs
ntr = length(dirdata.TrialNumber);
dirdata.postprob_ch_avg = nan(ntr,2);
sacc_time = [-50,100]/diff(t_mids(1:2));

for tr = 1:ntr
    
    % skip forced trials
    if bhvdata.trialtype(tr)==1
        continue
    end
    
    % sync to sacc1
    sacc1 = bhvdata.MLsacc(tr,1);
    if isnan(sacc1) | sacc1>1200
        continue;
    end
    
    [~, sacc1_zero] = min(abs(t_mids-sacc1));
    t_idx = sacc_time+sacc1_zero;
    dirdata.postprob_ch_avg(tr,1) = mean(dirdata.postprob_ch(tr,t_idx));
    
    % sync to sacc2
    sacc2 = bhvdata.MLsacc(tr,2);
    if isnan(sacc2) | sacc2>1200
        continue;
    end
    
    [~, sacc2_zero] = min(abs(t_mids-sacc2));
    t_idx = sacc_time+sacc2_zero;
    dirdata.postprob_ch_avg(tr,2) = mean(dirdata.postprob_ch(tr,t_idx));
    
end
dirdata.postprob_ch_avg = logit(dirdata.postprob_ch_avg);

figure; hold on
temp1 = dirdata.postprob_ch_avg(:,1);
temp1 = temp1(~isnan(temp1));
temp2 = dirdata.postprob_ch_avg(:,2);
temp2 = temp2(~isnan(temp2));
histogram(temp1,'Normalization','probability','BinWidth',.1)
histogram(temp2,'Normalization','probability','BinWidth',.1)

figure; hold on
temp3 = diff(dirdata.postprob_ch_avg,[],2);
temp3 = temp3(~isnan(temp3));
histogram(temp3,'Normalization','probability')

%% do anova on avg postprob ~ sacc# + saccloc
for S = 1:nsubj
    % restrict trials
    subject = subject_names{S};
    idx = find(strcmp(bhvdata.subject,subject_names{S}) & ...
        keep_tr);
    ntr = length(idx);
    
    % build anova vectors
    Y = [dirdata.postprob_ch_avg(idx,1);dirdata.postprob_ch_avg(idx,2)];
    saccID = repelem({'1st';'2nd'},ntr);
    saccloc = [bhvdata.saccloc_target(idx,1);bhvdata.saccloc_target(idx,2)];
    
    keep_rows = ~isnan(Y);
    Y = Y(keep_rows);
    saccID = saccID(keep_rows);
    saccloc = saccloc(keep_rows);
    
    Y_norm = (Y-mean(Y))./std(Y);
    
    [p,tbl,stats,terms] = anovan(Y_norm, {saccID,saccloc},...
        'model','interaction','varnames',{'saccID','saccloc'});
    figure; hold on
    saccID_opts = unique(saccID);
    saccloc_opts = [1,-1];
    track = nan(2,2,2);
    for i = 1:2
        for j = 1:2
            temp = Y_norm(strcmp(saccID,saccID_opts{i}) & saccloc==saccloc_opts(j));
            length(temp)
            track(i,j,:) = [mean(temp), std(temp)./sqrt(length(temp))];
        end
    end
    bar(track(:,:,1))
    errorbar([.85 1.15 1.85 2.15],[track(1,:,1),track(2,:,1)],[track(1,:,2),track(2,:,2)],'k.')
    set(gca,'XTick',1:2,'XTickLabel',saccID_opts)
    ylabel('norm logit posterior probability')
%     legend({'chosen','other'},'box','off','Location','southwest')
    title(subject)
end

%% do anova on delta avg postprob ~ type,
% where type = unch->ch > ch->ch > ch->unch
for S = 1:nsubj
    
    % restrict trials
    subject = subject_names{S};
    idx = find(strcmp(bhvdata.subject,subject_names{S}) & ...
        keep_tr);
    ntr = length(idx);
    
    % build anova vectors
    Y = diff(dirdata.postprob_ch_avg(idx,:),[],2);
    temp = bhvdata.saccloc_target(idx,:);
    temp(isnan(temp))=0;
    sacctype = mat2cell(num2str(temp),ones(ntr,1),5);
    sacctype_order = zeros(ntr,1);
    sacctype_order(ismember(sacctype,{'-1  1'})) = 1;
    sacctype_order(ismember(sacctype,{' 1  1'})) = 2;
   sacctype_order(ismember(sacctype,{' 1 -1'})) = 3;
    
    keep_rows = sacctype_order>0 & ~isnan(Y);
    
    Y = Y(keep_rows);
    sacctype_order = sacctype_order(keep_rows,:);
    
    sacctype = sacctype(keep_rows,:);
    
    Y_norm = (Y-mean(Y))./std(Y);
    
    [p,tbl,stats] = anova1(Y_norm,sacctype_order);
    
    figure; hold on
    track = nan(1,3,2);
    for i = 1:3
        temp = Y_norm(sacctype_order==i);
        track(1,i,:) = [mean(temp), std(temp)./sqrt(length(temp))];
    end
    bar(track(:,:,1))
    errorbar(1:3,track(:,:,1),track(:,:,2),'k.');
    set(gca,'XTick',1:3,'XTickLabel',...
        {'unch\rightarrowch','ch\rightarrowch','ch\rightarrowunch'});
    ylabel('norm logit \Delta posterior probability')
    title(subject)
  
end



