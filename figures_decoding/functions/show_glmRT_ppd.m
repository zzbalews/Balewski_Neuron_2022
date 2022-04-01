function show_glmRT_ppd(bhvdata,valdata,keep_tr)
%%
subject_names = unique(valdata.subject);
nsubj = length(subject_names);

t_mids = valdata.t_mids;
nmids = length(t_mids);

t_idx = find(~isnan(valdata.ch_ppd(1,:)));

thresh = 0.001;

pairs = [1,2;1,3;1,4;2,3;2,4;3,4];

figure;


for s = 1:nsubj
    
    % setup vars
vif = nan(2,nmids,7);
track_cpd = nan(2,nmids,7);
track_ifsig = nan(2,nmids,7);
track_B = nan(2,nmids,7);
track_ntr = nan(7,1);

  subject = subject_names{s};
    for i = 7%1:7
        
        % restrict trials
      
        if i<7 % split trials by pair
            tr_idx = find(strcmp(bhvdata.subject,subject) & ...
                keep_tr & ...
                min(bhvdata.valbin_expval,[],2)==pairs(i,1) & ...
                max(bhvdata.valbin_expval,[],2)==pairs(i,2));
        else % do all trials together
            tr_idx = find(strcmp(bhvdata.subject,subject) & ...
                keep_tr);
        end
        track_ntr(i)=length(tr_idx);
        
        % rt
        log10rt = log10(bhvdata.rt(tr_idx));
        Y_norm = (log10rt - mean(log10rt))./std(log10rt);
        
        % ppd
        ch = valdata.ch_ppd(tr_idx,:);
        unch = valdata.unch_ppd(tr_idx,:);
        
        % ppd na, randomly select option
        na = valdata.na_ppd(tr_idx,:,1);
        temp = randperm(length(tr_idx));
        temp = temp(1:round(length(tr_idx)/2));
        na(temp,:) = valdata.na_ppd(tr_idx(temp),:,2);
        
        % maxval
        maxval = max(bhvdata.valbin_expval(tr_idx,:),[],2);
        
        for t = t_idx
            
            % X = ch, unch
            temp = [ch(:,t),unch(:,t)];            
            X_norm = [ones(track_ntr(i),1),(temp-mean(temp))./std(temp)];
            
%                     % X = ch, na
%                     temp = [ch(:,t),na(:,t)];
%                     X_norm = [ones(track_ntr(i),1),(temp-mean(temp))./std(temp)];

%                     % X = unch, na
%                     temp = [unch(:,t),na(:,t)];%             
%                     X_norm = [ones(track_ntr(i),1),(temp-mean(temp))./std(temp)];
            
            [B_hat, ~, coeff_p, model_p, CPD] = run_glm_cpd(X_norm, Y_norm);
            track_cpd(:,t,i) = 100*CPD(2:3);
            track_ifsig(:,t,i) = coeff_p(2:3)<thresh & model_p<thresh;
            track_B(:,t,i) = B_hat(2:3);
            
            
        end
    end
    
    track_ifsig(track_ifsig==0) = NaN;
    
    clrs = [1 0 0 ; 0 0 1]; ll = {'ch','unch'};
%         clrs = [1 0 0; .5 .5 .5]; ll = {'ch','na'};
%         clrs = [0 0 1; .5 .5 .5]; ll = {'unch','na'};
    
    
    
    %% plot glm from all trials
%     figure;
%     subplot(2,1,1); hold on
    subplot(2,2,s+2); hold on
    p = plot(t_mids,track_cpd(:,:,7),'LineWidth',.5);
    for i = 1:2
        p(i).Color = [clrs(i,:) .3];
    end
    
    p = plot(t_mids,track_ifsig(:,:,7).*track_cpd(:,:,7),'LineWidth',2);
    for i = 1:2
        p(i).Color = clrs(i,:);
    end
    plot([0 0],[0 2],'k:')
    
    xlim([-500 1000])
    ylim([0 1.2])
    xlabel('time from pics (ms)')
    ylabel('% CPD')
    
    %     subplot(2,1,2); hold on
    subplot(2,2,s); hold on

    p = plot(t_mids,track_B(:,:,7),'LineWidth',.5);
    for i = 1:2
        p(i).Color = [clrs(i,:) .3];
    end
    
    p = plot(t_mids,track_ifsig(:,:,7).*track_B(:,:,7),'LineWidth',2);
    for i = 1:2
        p(i).Color = clrs(i,:);
    end
    plot([0 0],[-.3 .3],'k:')
    plot([-1000 2000],[0 0],'k:')
    
    %     legend(p,ll,'Location','best','box','off')
    
    xlim([-500 1000])
    ylim([-.15 .15])
    xlabel('time from pics (ms)')
    ylabel('\beta coefficient')
    
    text(-500, -.12,['  n=',num2str(track_ntr(7))])
    
%     
%     %% plot betas from trial groups
%     figure;
%     
%     for j = 1:6
%         
%         subplot(2,3,j); hold on
%         title([num2str(pairs(j,2)),' v ',num2str(pairs(j,1))])
%         p = plot(t_mids,track_B(:,:,j),'LineWidth',.5);
%         for i = 1:2
%             p(i).Color = [clrs(i,:) .3];
%         end
%         
%         p = plot(t_mids,track_ifsig(:,:,j).*track_B(:,:,j),'LineWidth',2);
%         for i = 1:2
%             p(i).Color = clrs(i,:);
%         end
%         plot([0 0],[-.3 .3],'k:')
%         plot([-1000 2000],[0 0],'k:')
%         
%         %     legend(p,ll,'Location','best','box','off')
%         
%         xlim([-500 1000])
%         ylim([-.3 .3])
%         xlabel('time from pics (ms)')
%         ylabel('\beta coefficient')
%         
%         text(-500, -.12,['  n=',num2str(track_ntr(j))])
%         
%     end
%     
%     
%     %% plot avg beta across trial groups
%     figure;
%     hold on
%     title(subject)
%     
%     for i = 1:2
%         shadedErrorBar(t_mids, mean(track_B(i,:,:),3),std(track_B(i,:,:),[],3)./sqrt(6),...
%             'lineprops',{'LineWidth',1,'Color',clrs(i,:)});
%     end
%     
%     plot([0 0],[-.3 .3],'k:')
%     plot([-1000 2000],[0 0],'k:')
%     xlim([-500 1000])
%     ylim([-.15 .15])
%     xlabel('time from pics (ms)')
%     ylabel('\beta coefficient')
    
end
end