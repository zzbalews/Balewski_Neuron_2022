function show_postprob_atsaccs(bhvdata,dirdata,keep_tr)



subject_names = unique(dirdata.subject);
nsubj = length(subject_names);

t_mids = dirdata.t_mids;
nmids = length(t_mids);

% relabel postprob by chosen direction
dirdata.postprob_ch = dirdata.postprob(:,:,1);
idxR = bhvdata.lever==1;
dirdata.postprob_ch(idxR,:) = dirdata.postprob(idxR,:,2);

% get maxvalue
bhvdata.maxvalbin = max(bhvdata.valbin_expval,[],2);
bhvdata.minvalbin = min(bhvdata.valbin_expval,[],2);
clrs_max = [0 0 0; 0 0 .75; 0 .2 1; .2 .5 1];
clrs_min = [0 0 0; .75 0 0; 1 .2 0; 1 .5 .2];

clrs_val = [0 0 0; 0 .7 0];
clrs_pattern = [.7 0 .7; 1 .5 0; .5 .5 .5; 0 .3 1];

% adjust trials by 1st saccade
ntr = size(bhvdata.TrialNumber,1);
dirdata.postprob_ch_sacc = nan(ntr,201,2);

t_mids_sacc = -500:5:500;

for tr = 1:size(bhvdata.TrialNumber,1)
    
    % skipp forced trials
    if bhvdata.trialtype(tr)==1
        continue
    end
    
    % sync to sacc1
    sacc1 = bhvdata.MLsacc(tr,1);
    if isnan(sacc1) | sacc1>1200
        continue;
    end
    
    [~, sacc1_zero] = min(abs(t_mids-sacc1));
    t_idx = [-100:100]+sacc1_zero;
    dirdata.postprob_ch_sacc(tr,:,1) = dirdata.postprob_ch(tr,t_idx);
    
    % sync to sacc2
    sacc2 = bhvdata.MLsacc(tr,2);
    if isnan(sacc2) | sacc2>1200
        continue;
    end
    
    [~, sacc2_zero] = min(abs(t_mids-sacc2));
    t_idx = [-100:100]+sacc2_zero;
    dirdata.postprob_ch_sacc(tr,:,2) = dirdata.postprob_ch(tr,t_idx);
    
end

%%

track_beta = nan(2,length(t_mids_sacc),2,2); % terms x time x which_sacc x subject
track_cpd = nan(size(track_beta));
track_ifsig = nan(size(track_beta));

track_unch_beta = nan(2,2); % terms x subject
track_unch_cpd = nan(2,2);
track_unch_p = nan(2,2);

dirdata.postprob_ch_sacc_resid = nan(size(dirdata.postprob_ch_sacc));

figure;
for S = 1:nsubj
    
    % restrict trials
    subject = subject_names{S};
    
    idx_subj = strcmp(bhvdata.subject,subject_names{S}) & ...
        keep_tr;
    
        idx_sacc1ch_sacc2ch = bhvdata.saccloc(:,1)==bhvdata.lever & ...
            bhvdata.saccloc(:,2)==bhvdata.lever;
        idx_sacc1ch_sacc2unch = bhvdata.saccloc(:,1)==bhvdata.lever & ...
            bhvdata.saccloc(:,2)==-bhvdata.lever;
        idx_sacc1ch_sacc2na = bhvdata.saccloc(:,1)==bhvdata.lever & ...
            isnan(bhvdata.saccloc(:,2));
        idx_sacc1unch = bhvdata.saccloc(:,1)==-bhvdata.lever;
    
        tr_groups = {idx_subj & idx_sacc1ch_sacc2ch, ...
            idx_subj & idx_sacc1ch_sacc2unch, ...
            idx_subj & idx_sacc1ch_sacc2na, ...
            idx_subj & idx_sacc1unch};
    
%     idx_val1 = bhvdata.maxvalbin==1;
%     idx_val2 = bhvdata.maxvalbin==2;
%     idx_val3 = bhvdata.maxvalbin==3;
%     idx_val4 = bhvdata.maxvalbin==4;
%     
%     tr_groups = {idx_subj & idx_val1, ...
%         idx_subj & idx_val2, ...
%         idx_subj & idx_val3, ...
%         idx_subj & idx_val4};
    
    Ngrp = length(tr_groups);
    
    med_sacc1 = nanmedian(bhvdata.MLsacc(idx_subj,1));
    med_sacc2 = nanmedian(bhvdata.MLsacc(idx_subj,2));
    
    % regressions & plot!
    for T = 1:2
               
        
        %% do regression! logit(post prob) ~ 1 + max value + saccade dir
        sacc_to_target = zeros(size(bhvdata.TrialNumber));
        sacc_to_target(bhvdata.saccloc(:,T)==bhvdata.lever) = 1;
        sacc_to_target(bhvdata.saccloc(:,T)==-bhvdata.lever) = -1;
        
        X = [bhvdata.maxvalbin, sacc_to_target];
        
        for t = 1:length(t_mids_sacc)
            
            idx_tr = idx_subj & ...
                ~isnan(bhvdata.saccloc(:,T)) & ...
                ~isnan(dirdata.postprob_ch_sacc(:,t,T));
            X_t = X(idx_tr,:);
            %             vif = max(diag(inv(corrcoef(X))))
            
            X_norm = [ones(sum(idx_tr),1), (X_t - mean(X_t))./std(X_t)];
            
            Y = logit(dirdata.postprob_ch_sacc(idx_tr,t,T));
            Y_norm = (Y - mean(Y))./std(Y);

            [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, Y_norm);
            
            track_beta(:,t,T,S) = B_hat(2:end);
            track_cpd(:,t,T,S) = CPD(2:end);
            track_ifsig(:,t,T,S) = coeff_p(2:end) < thresh & model_p < thresh;
            
        end
        
        
        %% do regression for residuals! logit(post prob) ~ 1 + max value --> residuals ~ saccade dir
        
        X = [bhvdata.maxvalbin];
        
        for t = 1:length(t_mids_sacc)
            
            idx_tr = idx_subj & ...
                ~isnan(bhvdata.saccloc(:,T)) & ...
                ~isnan(dirdata.postprob_ch_sacc(:,t,T));
            X_t = X(idx_tr,:);
                        
            X_norm = [ones(sum(idx_tr),1), (X_t - mean(X_t))./std(X_t)];
            
            Y = logit(dirdata.postprob_ch_sacc(idx_tr,t,T));
            Y_m = mean(Y);
            Y_sd = std(Y);
            Y_norm = (Y - Y_m)./Y_sd;
            
            [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, Y_norm);
            Y_fit = X_norm*B_hat;
            Y_resid_aftermaxval = Y_norm -Y_fit;
%             Y_fit = Y_sd*(X_norm*B_hat)+Y_m;
%             Y_resid_aftermaxval = invlogit(Y) - invlogit(Y_fit);
            
            dirdata.postprob_ch_sacc_resid(idx_tr,t,T) = Y_resid_aftermaxval;
            
        end
        
%         %% plot raw post prob
%         subplot(2,2,T + 2*(S-1)); hold on
%         
%         % ref lines
%         plot([0 0],[0 1],'k:')
%         plot([-2000 2000],[.5 .5],'k:')
%         
%         LL = [];
%         LL_label = {};
%         
%         for g = 1:Ngrp
%             
%             rowidx = tr_groups{g};
%             temp = dirdata.postprob_ch_sacc(rowidx,:,T);
%             a = shadedErrorBar(t_mids_sacc, nanmean(temp), ...
%                 nanstd(temp)./sqrt(sum(rowidx)),...
%                 'lineprops',{'Color',clrs_max(g,:)}); %clrs_pattern; clrs_max
%             
%             LL(g) = a.mainLine;
%             LL_label{g} = num2str(sum(rowidx));
%             
%             
%         end
%         
%         leg = legend(LL,LL_label,'box','off','Location','northwest');
%         title(leg,'n=');
%         xlim([-500 500])
%         ylim([.4 .8])
%         xlabel(['time from sacc #',num2str(T),' (ms)'])
%         ylabel({'direction decoder','chosen post. prob'})
%         title(subject)
%         
        
        %% plot logit(post prob) ~ 1 + max value --> residuals ~ saccade dir
        subplot(2,2,T + 2*(S-1)); hold on
        
        % ref lines
        
        plot([0 0],[-1 1],'k:')
        plot([-2000 2000],[0 0],'k:')
        
        LL = [];
        LL_label = {};
        
        for g = 1:Ngrp
            
            rowidx = tr_groups{g};
            temp = dirdata.postprob_ch_sacc_resid(rowidx,:,T);
            a = shadedErrorBar(t_mids_sacc, nanmean(temp), ...
                nanstd(temp)./sqrt(sum(rowidx)),...
                'lineprops',{'Color',clrs_pattern(g,:)}); %clrs_pattern; clrs_max
            
            LL(g) = a.mainLine;
            LL_label{g} = num2str(sum(rowidx));
            
            
        end
        
        leg = legend(LL,LL_label,'box','off','Location','northwest');
        title(leg,'n=');
        xlim([-500 500])
        ylim([-.3 .3])
        xlabel(['time from sacc #',num2str(T),' (ms)'])
        ylabel({'direction decoder','chosen post. prob'})
        title(subject)
        
        
        
    end 
    
    
    %% compare 100ms after 1st and 2nd saccade (if 1st=unchosen) to see if post prob increases
    t_idx = find(t_mids_sacc>0 & t_mids_sacc<=100);
    sacc1_100avg = nanmean(dirdata.postprob_ch_sacc(:,t_idx,1),2);
    sacc2_100avg = nanmean(dirdata.postprob_ch_sacc(:,t_idx,2),2);
    
    
    tr_idx = idx_subj & bhvdata.saccloc(:,1)==-bhvdata.lever & ...
        ~isnan(sacc1_100avg) & ~isnan(sacc2_100avg);
    

        Y = logit([sacc1_100avg(tr_idx);sacc2_100avg(tr_idx)]);
        Y_norm = (Y-mean(Y))./std(Y);
        X = [repmat(bhvdata.maxvalbin(tr_idx),2,1),repelem([-1;1],sum(tr_idx),1)];
        X_norm = [ones(2*sum(tr_idx),1), (X-mean(X))./std(X)];
        
         [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, Y_norm);
         
         track_unch_beta(:,S) = B_hat(2:end); 
track_unch_cpd(:,S) = CPD(2:end);
track_unch_p(:,S) = coeff_p(2:end);
model_p
        
end
track_ifsig = double(track_ifsig);
track_ifsig(track_ifsig==0) = NaN;

%% show regression
clrs_glm = [clrs_max(4,:);clrs_pattern(1,:)];

figure;
for S = 1:nsubj
    
    % restrict trials
    subject = subject_names{S};
    for T = 1:2
        subplot(2,2,T + 2*(S-1)); hold on
        plot([-500 500],[0 0],'k:')
        plot([0 0],[-1 1],'k:')
        LL = [];
       for x = 1:2
            plot(t_mids_sacc, track_beta(x,:,T,S), ...
                'Color', [.5 .5 .5],'LineWidth',1)
            a = plot(t_mids_sacc, track_beta(x,:,T,S).*track_ifsig(x,:,T,S),...
                'Color', clrs_glm(x,:), 'LineWidth',2);
            LL(x) = a;
        end  
        xlabel(['time from sacc #',num2str(T)])
        ylabel('\beta coefficient')
        title(subject)
        legend(LL,{'maxval','sacc to chosen'},'box','off','Location','southwest')
        xlim([-500 500])
        ylim([-.2 .2])
    end
end
    

%%

x = nan(2,6);

for S = 1:nsubj
    
    % restrict trials
    subject = subject_names{S};
    idx_subj = strcmp(bhvdata.subject,subject_names{S}) & ...
        keep_tr;
    
    idx_sacc1ch = bhvdata.saccloc(:,1)==bhvdata.lever;
    idx_sacc1unch = bhvdata.saccloc(:,1)==-bhvdata.lever;
    idx_sacc1na = isnan(bhvdata.saccloc(:,1));
    
    idx_sacc2ch = bhvdata.saccloc(:,2)==bhvdata.lever;
    idx_sacc2unch = bhvdata.saccloc(:,2)==-bhvdata.lever;
    idx_sacc2na = isnan(bhvdata.saccloc(:,2));
    
    
    x(S,:) = [sum(idx_subj & idx_sacc1ch & idx_sacc2ch),...
        sum(idx_subj & idx_sacc1ch & idx_sacc2unch),...
        sum(idx_subj & idx_sacc1ch & idx_sacc2na),...
        sum(idx_subj & idx_sacc1unch & idx_sacc2ch),...
        sum(idx_subj & idx_sacc1unch & idx_sacc2unch),...
        sum(idx_subj & idx_sacc1na)]./sum(idx_subj);
    
    round(x,2)
    
    
end
figure;
b = bar(x(:,end:-1:1),'stacked');
b(6).FaceColor = clrs_pattern(1,:);
b(5).FaceColor = clrs_pattern(2,:);
b(4).FaceColor = clrs_pattern(3,:);
b(3).FaceColor = clrs_pattern(4,:);
b(2).FaceColor = [0 0 0];
b(1).FaceColor = [0 0 0];
for j = 1:6
    b(j).FaceAlpha = .7;
end

end

