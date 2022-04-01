function show_responsetime_vsvalue(bhvdata,variable,which_value)
%show 'variable' as a function of 'which_value', for only 'trialtype' trials



%% show subject results (separte plots)
figure('Name',[variable]);

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
    
    switch which_value
        case 'maxval'
            val = max(bhvdata.subjval_expval,[],2);
            add_forced = 1;
        case 'valdiff'
            bhvdata.subjval_expval(isnan(bhvdata.subjval_expval)) = 0;
            val = abs(diff(bhvdata.subjval_expval,[],2));
            add_forced = 0;
    end
    
    edges = linspace(min(val),max(val),6);
    midpoints = mean([edges(1:end-1); edges(2:end)]);
    
    val_disc = midpoints(discretize(val,edges))';
    
    switch variable
        case 'RT'
            y = bhvdata.rt;
        case 'sacctime'
            y = bhvdata.MLsacc;
    end
    
    track_variable_free = nan(length(midpoints),2);
        track_variable_forced = nan(length(midpoints),2);

    for i = 1:length(midpoints)
        % free trials
        idx = val_disc==midpoints(i) & tr_subj_free;
        N = sum(~isnan(y(idx)));
        track_variable_free(i,:) = [nanmean(y(idx)), nanstd(y(idx))./sqrt(N)];
        
        if add_forced
        % forced trials
        idx = val_disc==midpoints(i) & tr_subj_forced;
        N = sum(~isnan(y(idx)));
        track_variable_forced(i,:) = [nanmean(y(idx)), nanstd(y(idx))./sqrt(N)];
        end
        
    end
    
    % plot
    subplot(1,2,s); hold on
    a=shadedErrorBar(midpoints,track_variable_free(:,1),track_variable_free(:,2),...
        'lineprops',{});
    if add_forced
    b=shadedErrorBar(midpoints,track_variable_forced(:,1),track_variable_forced(:,2),...
        'lineprops',{});
    end
    xlabel(which_value)
    ylabel(variable)
     
    title(subject)
    if add_forced
    legend([a.mainLine,b.mainLine],{'free','forced'},'box','off')
    end
    
    
%     %% check glm
%     valmax = max(bhvdata.subjval_expval,[],2);
%     valdiff = abs(diff(bhvdata.subjval_expval,[],2));
%     
%     X = [valmax,valdiff];
%     X = X(tr_subj_free,:);
%     X_norm = [ones(size(X,1),1),(X-mean(X))./std(X)];
%     
%     Y = bhvdata.rt(tr_subj_free);
%     Y_norm = (Y-mean(Y))./std(Y);
%     
%     [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, Y_norm);
%     
end





end