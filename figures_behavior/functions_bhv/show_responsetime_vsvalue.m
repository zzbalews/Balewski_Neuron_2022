function show_responsetime_vsvalue(bhvdata,variable,which_value)
%show 'variable' as a function of 'which_value', for only 'trialtype' trials

switch which_value
    case 'subjval'
        bhvdata.VALUE = bhvdata.subjval_expval;
    case 'valbin'
        bhvdata.VALUE = bhvdata.valbin_expval;
end


%% show subject results (separte plots)
figure('Name',[variable]);

subject_names = unique(bhvdata.subject);



for s = 1:length(subject_names)
    
    % restrict to ('trialtype') trials from subject
    subject = subject_names{s};
    
    switch variable
        case 'RT'
            tr_subj_free = strcmp(bhvdata.subject,subject) & ...
                bhvdata.trialtype==2 & ...
                bhvdata.lever~=0;
        case 'sacctime'
            tr_subj_free = strcmp(bhvdata.subject,subject) & ...
                bhvdata.trialtype==2 & ...
                bhvdata.lever~=0 & ...
                ~isnan(bhvdata.saccloc(:,1));
    end
    
    valmax = max(bhvdata.VALUE(tr_subj_free,:),[],2);
%     valmax_m = mean(valmax);
%     valmax_sd = std(valmax);
%     valmax_norm = (valmax-valmax_m)./valmax_sd;
    
    valdiff = abs(diff(bhvdata.VALUE(tr_subj_free,:),[],2));
%     valdiff_m = mean(valdiff);
%     valdiff_sd = std(valdiff);
%     valdiff_norm = (valdiff-valdiff_m)./valdiff_sd;
    
%     vif = max(diag(inv(corrcoef([valmax_norm,valdiff_norm]))));
    
    switch variable
        case 'RT'
%             y = log10(bhvdata.rt(tr_subj_free));
            y = bhvdata.rt(tr_subj_free)/1000;
        case 'sacctime'
%             y = log10(bhvdata.MLsacc(tr_subj_free,1));
            y = bhvdata.MLsacc(tr_subj_free,1)/1000;
    end
    
%     y_m = nanmean(y);
%     y_sd = nanstd(y);
%     y_norm = (y-y_m)./y_sd;
    
%     % full model
%     X_norm = [ones(sum(tr_subj_free),1),valmax_norm,valdiff_norm];
%     [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, y_norm);
%     disp(subject)
%     disp(['... GLM: log10(',num2str(variable),') ~ 1 + maxval + valdiff'])
%     disp(['... maxval beta = ',num2str(round(B_hat(2),2)),...
%         ', p=',num2str(coeff_p(2)),...
%         ', CPD=', num2str(round(100*CPD(2),1))])
%     disp(['... valdiff beta = ',num2str(round(B_hat(3),2)),...
%         ', p=',num2str(coeff_p(3)),...
%         ', CPD=', num2str(round(100*CPD(3),1))])
    
    % y ~ valmax
%     X_norm = [ones(sum(tr_subj_free),1),valmax_norm];
%     [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, y_norm);
%     y_fit = y_sd*(X_norm*B_hat)+y_m;
%     y_resid_aftervalmax = 10.^y - 10.^y_fit;
%     
%     % y ~ valdiff
%     X_norm = [ones(sum(tr_subj_free),1),valdiff_norm];
%     [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, y_norm);
%     y_fit = y_sd*(X_norm*B_hat)+y_m;
%     y_resid_aftervaldiff = 10.^y - 10.^y_fit;
    
    
    % discretize x axes for nice plots
    nmids = 4;
    edges = linspace(min(valmax),max(valmax),nmids+1);
    valmax_midpoints = mean([edges(1:end-1); edges(2:end)]);
    valmax_disc = valmax_midpoints(discretize(valmax,edges))';
    
    edges = linspace(min(valdiff),max(valdiff),nmids+1);
    valdiff_midpoints = mean([edges(1:end-1); edges(2:end)]);
    valdiff_disc = valdiff_midpoints(discretize(valdiff,edges))';
    
    % track avg residuals
    track_y_resid_aftervalmax = nan(nmids,2);
    track_y_resid_aftervaldiff = nan(nmids,2);
    
    for i = 1:nmids
        
        % y ~ valmax --> residuals ~ valdiff
        idx = valdiff_disc==valdiff_midpoints(i);
        N = sum(~isnan(y(idx)));
        
        track_y_resid_aftervalmax(i,:) = [...
            nanmean(y(idx)),...
            nanstd(y(idx))./sqrt(N)];
        
        % y ~ valdiff --> residuals ~ valmax
        idx = valmax_disc==valmax_midpoints(i);
        N = sum(~isnan(y(idx)));
        
        track_y_resid_aftervaldiff(i,:) = [...
            nanmean(y(idx)),...
            nanstd(y(idx))./sqrt(N)];
    end
    
    track_y = nan(0,2);
    if strcmp(variable,'RT')
        style = '-';
    elseif strcmp(variable,'sacctime')
        style = '-'; %':';
    end
    
    % plot: y ~ valmax --> residuals ~ valdiff
    ax1 = subplot(2,2,s); hold on
    a=shadedErrorBar(valdiff_midpoints,track_y_resid_aftervalmax(:,1),...
        track_y_resid_aftervalmax(:,2),...
        'lineprops',{'Color',[0 .8 .4],'LineWidth',2,'LineStyle',style});
    xlabel(['\Delta',which_value])
    ylabel({[],[variable,' (s)']})
    title(subject)
    temp = valdiff_midpoints([1,end]);
    xlim(.2*diff(temp)*[-1 1]  + temp)
    temp = get(gca,'YLim');
    track_y = cat(1,track_y,temp);
    set(gca,'XTick',0:3)
    
    
    % plot: y ~ valdiff --> residuals ~ valmax
    ax2 = subplot(2,2,2+s); hold on
    a=shadedErrorBar(valmax_midpoints,track_y_resid_aftervaldiff(:,1),...
        track_y_resid_aftervaldiff(:,2),...
        'lineprops',{'Color',[0 .8 .4],'LineWidth',2,'LineStyle',style});
    xlabel(['max ',which_value])
    ylabel({[],[variable,' (s)']})
    title(subject)
    temp = valmax_midpoints([1,end]);
    xlim(.2*diff(temp)*[-1 1]  + temp)
    temp = get(gca,'YLim');
    track_y = cat(1,track_y,temp);
    set(gca,'XTick',1:4)
    
    Y = [50*floor(min(1000*track_y(:,1))/50),50*ceil(max(1000*track_y(:,2))/50)]/1000;
    for i = [0,2]
        subplot(2,2,i+s)
        ylim(Y);
    end
    
    
    
end



end