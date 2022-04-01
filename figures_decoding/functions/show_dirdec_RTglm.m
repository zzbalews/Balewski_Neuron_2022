function show_dirdec_RTglm(bhvdata,dirdata,keep_tr,thresh)

subject_names = unique(dirdata.subject);
nsubj = length(subject_names);

t_mids = dirdata.t_mids;
nmids = length(t_mids);


% relabel by chosen direction
dirdata.postprob_ch = dirdata.postprob(:,:,1);
idxR = bhvdata.lever==1;
dirdata.postprob_ch(idxR,:) = dirdata.postprob(idxR,:,2);


clrs = [1 .6 0];


pairs = [1,2;1,3;1,4;2,3;2,4;3,4];
figure;

for S = 1:nsubj
    
    % setup vars
    vif = nan(1,nmids,7);
    track_cpd = nan(1,nmids,7);
    track_ifsig = nan(1,nmids,7);
    track_B = nan(1,nmids,7);
    track_ntr = nan(7,1);
    track_pval = nan(1,nmids,7);
    subject = subject_names{S};
    for i = 7%1:7
        
        
        % restrict trials
        if i<7 % split trials by pair
            idx = find(strcmp(bhvdata.subject,subject) & ...
                keep_tr & ...
                min(bhvdata.valbin_expval,[],2)==pairs(i,1) & ...
                max(bhvdata.valbin_expval,[],2)==pairs(i,2));
        else % do all trials together
            idx = find(strcmp(bhvdata.subject,subject) & ...
                keep_tr);
        end
        track_ntr(i)=length(idx);
        
        %     % relabel by chosen direction
        %     chosen_dirdec = dirdata.postprob(idx,:,1);
        %     idxR = bhvdata.lever(idx)==1;
        %     chosen_dirdec(idxR,:) = dirdata.postprob(idx(idxR),:,2);
        
        % do GLM: RT ~ postprob
        
        maxval = max(bhvdata.valbin_expval(idx,:),[],2);
        maxval_norm = (maxval - mean(maxval))./std(maxval);
        
        Y = log10(bhvdata.rt(idx));
        Y_norm = (Y - mean(Y))./std(Y);
        
        t_idx = find(t_mids>=-500 & t_mids<=1500);
        
        
        for t = t_idx
            X = logit(dirdata.postprob_ch(idx,t));
            
            if i<7
                X_norm = [ones(track_ntr(i),1),(X - mean(X))./std(X)];
            else
                X_norm = [ones(track_ntr(i),1),(X - mean(X))./std(X), maxval_norm];
            end
            
            %         vif(S,t) = max(diag(inv(corrcoef(X_norm(:,2:3)))));
            
            
            [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, Y_norm);
            track_B(1,t,i) = B_hat(2);
            track_cpd(1,t,i) = CPD(2);
            track_ifsig(1,t,i) = coeff_p(2)<thresh & model_p<thresh;
            track_pval(1,t,i) = coeff_p(2);
            
        end
    end
    
    
    track_ifsig(track_ifsig==0) = NaN;
    
    
    %% show regression results (all trials)
    
%     figure;
    
[~,peak] = min(track_B(1,:,7));
        disp([subject,': RT ~ 1 + postprob + maxvalue'])
        disp(['    @peak=',num2str(1+t_mids(peak)),...
            ', B=',num2str(track_B(1,peak,7)),...
            ', %CPD=',num2str(100*track_cpd(1,peak,7)),...
            ', p=',num2str(track_pval(1,peak,7))])

%     subplot(2,1,1); hold on
    subplot(2,2,S); hold on
    
    plot(t_mids,track_B(:,:,7),'Color',[.5 .5 .5],'LineWidth',1);
    plot(t_mids,track_ifsig(:,:,7).*track_B(:,:,7),'Color',clrs,'LineWidth',2);
    set(gca,'XTick',-500:500:1500,'XTickLabel',-.5:.5:1.5,...
        'YTick',[-.1 0])
    
    xlim([-500 1500])
    ylim([-.175,.05])
    
    plot([0 0 ],[-1 1],'k:')
    plot([-2000 2000],[0 0],'k:')
    
    
    xlabel('time from pics (ms)')
    ylabel('\beta coefficient')
    
%     subplot(2,1,2); hold on
    subplot(2,2,2+S); hold on
    
    plot(t_mids,100*track_cpd(:,:,7),'Color',[.5 .5 .5],'LineWidth',1);
    plot(t_mids,100*track_ifsig(:,:,7).*track_cpd(:,:,7),'Color',clrs,'LineWidth',2);
    set(gca,'XTick',-500:500:1500,'XTickLabel',-.5:.5:1.5,...
        'YTick',0:2)
    xlim([-500 1500])
    ylim([0 2.75])
    
    plot([0 0 ],[0 100],'k:')
    
    xlabel('time from pics (ms)')
    ylabel('%CPD')
    
    
%     %% plot betas from trial groups
%     figure;
%     for j = 1:6
%         subplot(2,3,j); hold on
%         title([num2str(pairs(j,2)),' v ',num2str(pairs(j,1))])
%         p = plot(t_mids,track_B(:,:,j),'LineWidth',.5,'Color',[clrs .3]);
%         p = plot(t_mids,track_ifsig(:,:,j).*track_B(:,:,j),'LineWidth',2,'Color',clrs);
%         plot([0 0],[-.3 .3],'k:')
%         plot([-1000 2000],[0 0],'k:')
%         
%         xlim([-500 1000])
%         ylim([-.3 .3])
%         xlabel('time from pics (ms)')
%         ylabel('\beta coefficient')
%         text(-500, -.12,['  n=',num2str(track_ntr(j))])
%         
%     end
%     
%     %% plot avg beta across trial groups
%     figure; hold on
% %     title(subject)
%     shadedErrorBar(t_mids, mean(track_B(1,:,:),3), std(track_B(1,:,:),[],3)./sqrt(6),...
%         'lineprops',{'LineWidth',1,'Color',clrs});
%     plot([0 0],[-.3 .3],'k:')
%     plot([-1000 2000],[0 0],'k:')
%     xlim([-500 1000])
%     ylim([-.2 .2])
%     xlabel('time from pics (ms)')
%     ylabel('\beta coefficient')
%     
    
end

