function show_dirdec_firststate(bhvdata,dirdata,keep_tr)

%%
subject_names = unique(bhvdata.subject);
nsubj = length(subject_names);

t_mids = dirdata.t_mids;
[~,t0] = min(abs(t_mids-0));
figure;
for S = 1:nsubj
    % restrict trials
    subject = subject_names{S};
    idx = find(strcmp(bhvdata.subject,subject) & ...
        keep_tr);
    
    ntr = length(idx);
    
    % get 1st saccade
    sacc = bhvdata.MLsacc(idx,1);
    
    % get onset of first direction "state"
    state = dirdata.postprob_ch_disc(idx,t0:end);
    onset = nan(ntr,1);
    cn_dir = nan(ntr,1);
    for tr = 1:ntr
        temp = find(state(tr,:)~=0,1);
        if ~isempty(temp)
            onset(tr) = t_mids(t0-1+temp) ;
            cn_dir(tr) = state(tr,temp);
        end
    end
    
    %% hist 1st state onset & 1st sacc onset
    subplot(1,2,S); hold on
    histogram(onset,'Normalization','probability','BinWidth',10,...
        'FaceColor',[0 .4 .8],'EdgeAlpha',0)
    
    histogram(sacc,'Normalization','probability','BinWidth',10,...
        'DisplayStyle','stairs','EdgeColor',[0 .4 .8])
    disp(['state: ',num2str(nanmedian(onset))])
    disp(['sacc: ',num2str(nanmedian(sacc))])
    
    xlim([0 600])
    xlabel('Time from pictures (s)')
    ylabel('norm count')
    set(gca,'XTick',0:300:2000,'XTickLabel',0:.3:2)
    [H,P,CI,STATS] = ttest(onset,sacc);
    disp(['state vs sacc onset: t(',num2str(STATS.df),')=',num2str(STATS.tstat),',p=',num2str(P)])
    

    keep = ~isnan(onset) & ~isnan(sacc);
    [rho,p] = corr(onset(keep),sacc(keep));
    disp(['state vs sacc onset: \rho = ',num2str(rho),', p = ',num2str(p)])
    
    %% compare direction for 1st saccade and 1st state
    saccdir = 2*(bhvdata.saccloc(idx,1)==bhvdata.lever(idx))-1;
    saccdir(isnan(bhvdata.saccloc(idx,1))) = 0;
    
    temp = [saccdir,cn_dir];
    temp(isnan(temp)) = 0;
    
    disp([subject,': state/sacc direction agreement = ',num2str(mean(temp(:,1)==temp(:,2)))])
    
    groups = {[1 1],[1 -1];[-1 1],[-1 -1]};
    observed = nan(2,2);
    for i = 1:2
       for j = 1:2
           tr = sum(temp==groups{i,j},2)==2;
            observed(i,j) = sum(tr);
       end
    end
    
    [chi2,df,chi2p] = do_chi2(observed);
    disp(['   chi2 test: chi2 = ',num2str(chi2),', df = ',num2str(df), 'p = ',num2str(chi2p)])
    
    %%
   
%     legend({'onset dir state','1st saccade'},'box','off')
%     xlim([0 600])
%     ylim([0 .22])
%     
%     set(gca,'XTick',0:300:600,'XTickLabel',0:.3:.6,'YTick',0:.1:.2)
%     
%     xlabel('time from pics on (s)')
%     ylabel('norm. trial count')
    
%     %% scatter 1st state onset vs 1st sacc onset
%     subplot(1,2,S);
%     scatter(onset,sacc,'filled','MarkerFaceAlpha',.1)
%     xlabel('onset')
%     ylabel('sacc')
%     
%     maxval = max(bhvdata.valbin_expval(idx,:),[],2);
%     keep_rows = ~isnan(onset) & ~isnan(sacc) & ~isnan(maxval);
%     
%     Y = sacc(keep_rows);
%     Y_norm = (Y-mean(Y))./std(Y);
%     
%     X = [onset(keep_rows), maxval(keep_rows)];
%     X_norm = [ones(sum(keep_rows),1), (X-mean(X))./std(X)];
%     
%     [ B_hat, B_se, coeff_p, model_p, CPD] = run_glm_cpd(X_norm, Y_norm);
%     disp([subject,' GLM: 1st sacc ~ state onset'])
%     disp(['   onset: p=',num2str(coeff_p(2)),', B=',num2str(round(B_hat(2),2)),', CPD=',num2str(round(100*CPD(2),2))])
%      
    
%     %% GLM: 1st state onset ~ maxval + log10(RT)
%     
%     % get terms & kick out NaN rows
%     Y = onset;
%     maxval = max(bhvdata.valbin_expval(idx,:),[],2);
%     log10rt = log10(bhvdata.rt(idx));
%     
%     keep_rows = ~isnan(onset) & ~isnan(maxval) & log10rt>2 & ~isnan(sacc);
%     
%     Y = Y(keep_rows);
%     maxval = maxval(keep_rows);
%     log10rt = log10rt(keep_rows);
%     sacc = sacc(keep_rows);
%     
%     Y_m = mean(Y);
%     Y_sd = std(Y);
%     Y_norm = (Y - Y_m)./Y_sd;
%     
%     
%     % regression with both terms
% %     X = maxval;
%     X = [maxval, log10rt,sacc];
%        vif = max(diag(inv(corrcoef(X))))
%     X_norm = [ones(sum(keep_rows),1), (X - mean(X))./std(X)];
%     
%     [ B_hat, B_se, coeff_p, model_p, CPD] = run_glm_cpd(X_norm, Y_norm);
% %     disp([subject,' GLM: state onset ~ maxval + log10(RT)'])
%     disp([subject,' GLM: state onset ~ maxval'])
% disp(['   maxval: p=',num2str(coeff_p(2)),', B=',num2str(round(B_hat(2),2)),', CPD=',num2str(round(100*CPD(2),2))])
%     disp(['   log10(RT): p=',num2str(coeff_p(3)),', B=',num2str(round(B_hat(3),2)),', CPD=',num2str(round(100*CPD(3),2))])
%     disp(['   sacc: p=',num2str(coeff_p(4)),', B=',num2str(round(B_hat(4),2)),', CPD=',num2str(round(100*CPD(4),2))])
%     
%     
%     track_Y = nan(4,2);
%     
%     for i = 1:4
%         
%         temp = maxval==i;
%         N = sum(temp);  
%         
%         track_Y(i,:) = [...
%             nanmean(Y(temp)), ...
%             nanstd(Y(temp))./sqrt(N)];
%         
%     end
%     
%     subplot(1,2,S); hold on
%     a = shadedErrorBar(1:4,...
%         track_Y(:,1),...
%         track_Y(:,2),...
%         'lineprops',{'Color',[0 .4 .8]});
%     set(gca,'YTick',150:50:350,'YTickLabel',.150:.05:.35)
%     ylim([180 320])
%     xlim([.5 4.5])
%     xlabel('maxval')
%     ylabel('1st direction state onset (s)')
    
end


end

