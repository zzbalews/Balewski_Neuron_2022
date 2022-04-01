function show_RT_glm(bhvdata,valdata,dirdata,keep_tr,thresh)
%%
subject_names = unique(dirdata.subject);
nsubj = length(subject_names);

t_mids = dirdata.t_mids;
nmids = length(t_mids);

% relabel dir dec by chosen direction
dirdata.postprob_ch = dirdata.postprob(:,:,1);
idxR = bhvdata.lever==1;
dirdata.postprob_ch (idxR,:) = dirdata.postprob(idxR,:,2);

nterms = 3;
clrs = [.2 .8 .2;%1 .6 0;
    1 0 0;
    0 0 1];
ll = {'CNdir','OFCval ch','OFCval unch'};

pairs = [1,2;1,3;1,4;2,3;2,4;3,4];

for S = 1:nsubj
    
    % setup vars
    vif = nan(1,nmids,7);
    track_cpd = nan(nterms,nmids,7);
    track_ifsig = nan(nterms,nmids,7);
    track_B = nan(nterms,nmids,7);
    track_ntr = nan(7,1);
    
    subject = subject_names{S};
    
    for i = 7 %1:7
        
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
        
        track_ntr(i) = length(idx);
        
        
        
        % max val
        maxval = max(bhvdata.valbin_expval(idx,:),[],2);
        
        % RT
        Y = log10(bhvdata.rt(idx));
        Y_norm = (Y - mean(Y))./std(Y);
        
        % do GLM over time
        t_idx = find(t_mids>=-500 & t_mids<=1500);
        
        for t = t_idx
            
            % CN direction (ch)
            x1 = logit(dirdata.postprob_ch(idx,t));
            
            % OFC value (ch)
            x2 = valdata.ch_ppd(idx,t);
            
            % OFC value(unch)
            x3 = valdata.unch_ppd(idx,t);
            
            % norm X
            X = [x1,x2,x3];
%             X = [x1,x2];
            vif(1,t,i) = max(diag(inv(corrcoef(X))));
            X_norm = [ones(track_ntr(i),1), (X-mean(X))./std(X)];
            
            [ B_hat, ~, coeff_p, model_p, CPD, ~, ~] = run_glm_cpd(X_norm, Y_norm);
            track_B(:,t,i) = B_hat(2:end);
            track_cpd(:,t,i) = 100*CPD(2:end);
            track_ifsig(:,t,i) = coeff_p(2:end)<thresh & model_p<thresh;
            
        end
    end
    
    track_ifsig(track_ifsig==0) = NaN;
    
    
    %% plot glm from all trials
    figure;
    subplot(2,1,1); hold on
    
     plot([-2000 2000],[0 0],'k:')
    plot([0 0],[-.3 .3],'k:')
    p = plot(t_mids,track_B(:,:,7),'Color',[.5 .5 .5 .5],'LineWidth',1);
    for i = 1:nterms
        p(i).Color = [clrs(i,:),.5];
    end
    p = plot(t_mids,track_ifsig(:,:,7).*track_B(:,:,7),'LineWidth',2);
    for i = 1:nterms
        p(i).Color = clrs(i,:);
    end
    xlim([-500 1000]);
    ylim([-.3 .3])
    xlabel('time from pics (ms)')
    ylabel('\beta coefficients')
    %             legend(p,ll(1:nterms),'box','off','Location','best')
    
    temp = track_cpd(:,:,7).*track_ifsig(:,:,7);
   temp = cat(1,nan(1,nmids),temp);
    [~,row] = max(temp);
best_term = row-1;


for i = 1:nterms
    temp = double(best_term==i);
    temp(temp==0) = NaN;
plot(t_mids,temp*.17,'Color',clrs(i,:),'LineWidth',3);
end
    
    subplot(2,1,2); hold on
    plot([0 0],[0 7],'k:')
    p = plot(t_mids,track_cpd(:,:,7),'Color',[.5 .5 .5 .5],'LineWidth',1);
    p = plot(t_mids, track_ifsig(:,:,7).*track_cpd(:,:,7),'LineWidth',2);
    for i = 1:nterms
        p(i).Color = clrs(i,:);
    end
    xlim([-500 1000]);
    ylim([0 7])
    xlabel('time from pics (ms)')
    ylabel('%CPD')
    
   
% % %     
% % %     %% plot betas from trial groups
% % %     figure; 
% % %     for j = 1:6
% % %         
% % %         subplot(2,3,j); hold on
% % %         
% % %                 title([num2str(pairs(j,2)),' v ',num2str(pairs(j,1))])
% % % 
% % %                 
% % %                 p = plot(t_mids,track_B(:,:,j),'LineWidth',.5);
% % %         for i = 1:3
% % %             p(i).Color = [clrs(i,:) .3];
% % %         end
% % %         
% % %          p = plot(t_mids,track_ifsig(:,:,j).*track_B(:,:,j),'LineWidth',2);
% % %         for i = 1:3
% % %             p(i).Color = clrs(i,:);
% % %         end
% % %         plot([0 0],[-.3 .3],'k:')
% % %         plot([-1000 2000],[0 0],'k:')
% % %         xlim([-500 1000])
% % %         ylim([-.3 .3])
% % %         xlabel('time from pics (ms)')
% % %         ylabel('\beta coefficient')
% % %         
% % %           text(-500, -.12,['  n=',num2str(track_ntr(j))])
% % %         
% % %     end
% % % 
% % % %% plot avg beta across trial groups
% % % 
% % % figure; hold on
% % % title(subject)
% % % 
% % % for i = 1:3
% % %     shadedErrorBar(t_mids, mean(track_B(i,:,:),3),...
% % %         std(track_B(i,:,:),[],3)./sqrt(6),...
% % %         'lineprops',{'LineWidth',1,'Color',clrs(i,:)});
% % % end
% % % 
% % %   plot([0 0],[-.3 .3],'k:')
% % %     plot([-1000 2000],[0 0],'k:')
% % %     xlim([-500 1000])
% % %     ylim([-.15 .15])
% % %     xlabel('time from pics (ms)')
% % %     ylabel('\beta coefficient')
% % %     
end
