function show_dirdec_bymaxval(bhvdata,dirdata,keep_tr,thresh,add_eyes)

subject_names = unique(dirdata.subject);
nsubj = length(subject_names);

t_mids = dirdata.t_mids;
nmids = length(t_mids);

clrs = [linspace(.5,0,4)',linspace(0,.7,4)',linspace(.5,1,4)'];

tr_match1 = bhvdata.saccloc(:,1)==bhvdata.lever;
tr_opp1 = bhvdata.saccloc(:,1)==-bhvdata.lever;

tr_match2 = bhvdata.saccloc(:,2)==bhvdata.lever & tr_match1;
tr_opp2 = bhvdata.saccloc(:,2)==-bhvdata.lever & tr_match1;

clrs_glm = [clrs(4,:,:);1 .7 .3];

if add_eyes
    nterms = 2;
    vif = nan(2,nmids);
    t_eyes = -2000:2000;
    juice = [1075,675];
    glmvarnames = {'maxvalue','gaze location'};
    ROWS = 2;
else
    nterms = 1;
    glmvarnames = {'maxvalue'};
    ROWS = 1;
end



beta = nan(2,nmids,nterms);
cpd = nan(2,nmids,nterms);
ifsig = nan(2,nmids,nterms);
pval = nan(2,nmids,nterms);

figure;
for S = 1:nsubj
    
      % restrict trials
    subject = subject_names{S};
    idx = find(strcmp(bhvdata.subject,subject_names{S}) & ...
        keep_tr);
    
    ntr = length(idx);
    
    % relabel by chosen direction
    chosen_dirdec = dirdata.postprob(idx,:,1);
    idxR = bhvdata.lever(idx)==1;
    chosen_dirdec(idxR,:) = dirdata.postprob(idx(idxR),:,2);
    
    % split by maxval levels
    maxval = max(bhvdata.valbin_expval(idx,:),[],2);
    
%     % rt
%     rt = bhvdata.rt(idx);
%     
    % plot!
    if add_eyes
    tr_groups = {tr_match1(idx),tr_opp1(idx)};
    else
        tr_groups = {ones(size(idx))};
    end
    %         tr_groups = {tr_match2(idx),tr_opp2(idx)};
    sacc1time = bhvdata.MLsacc(idx,1);
    sacc2time = bhvdata.MLsacc(idx,2);
    
    for j = 1:ROWS
        subplot(ROWS,2,S+(j-1)*2); hold on
        plot([0 0],[0 1],'k:')
        plot([-2000 2000],[.5 .5],'k:')
        
%         % sacc time
%         temp = median(sacc1time(tr_groups{j}==1));
%         plot(temp*[1 1],[0 1],'r')
%         temp = nanmedian(sacc2time(tr_groups{j}==1));
%         plot(temp*[1 1],[0 1],'r')
        
%         % rt
%         
%         plot(median(rt(tr_groups{j}))*[1 1],[0 1],'g')
        ll = [];
        for v = 1:4
            rowidx = maxval==v & tr_groups{j};
            aa = shadedErrorBar(t_mids,mean(chosen_dirdec(rowidx,:)),...
                std(chosen_dirdec(rowidx,:))/sqrt(sum(rowidx)),...
                'lineprops',{'Color',clrs(v,:),'LineWidth',2});
            ll(v) = aa.mainLine;
            
            %         plot(median(rt(rowidx))*[1 1],[0 1],'Color',[0 v 0]/4)
            
            
        end
        
        
        
        
        % format
        xlim([-500 1500])
        ylim([.4 .8])
        set(gca,'XTick',-500:500:1500,'XTickLabel',-.5:.5:1.5)
        xlabel('time from pics (s)')
        ylabel({'direction decoder','chosen posterior probability'})
        
        leg = legend(ll(end:-1:1),{'4','3','2','1'},'box','off','Location','northeast');
        title(leg,'max value')
    end
    
    
    
    %% do GLM: post prob ~ maxval + gaze loc
    
    
    if add_eyes
        rt = bhvdata.rt(idx);
        
        % if looking at chosen pic?
        gaze = bhvdata.location_target(idx,:);
        ll = bhvdata.lever(idx)==-1;
        gaze(ll,:) = -gaze(ll,:);
        for tr = 1:size(gaze,1)
            bound = min(round(rt(tr)+2000+juice(S)),4000);
            gaze(tr,bound:end) = NaN;
        end
        t_eyes_downsample = ismember(t_eyes,t_mids);
        gaze = gaze(:,t_eyes_downsample);
        
    end
    
    
    t_idx = find(t_mids>=-500 & t_mids<=1500);
    
    for t = t_idx
        
        if add_eyes
            tr_in_slice = find(~isnan(gaze(:,t)));
            eyes_in_slice = length(unique(gaze(tr_in_slice,t)))>1;
        else
            tr_in_slice = 1:length(maxval);
            eyes_in_slice = 0;
        end
        
        ntr = length(tr_in_slice);
        
        if ntr==0
            continue
        end
        
        if eyes_in_slice
            X = [maxval(tr_in_slice),gaze(tr_in_slice,t)];
            vif(S,t) = max(diag(inv(corrcoef(X))));
        else
            X = [maxval(tr_in_slice)];
        end
        nX = size(X,2);
        X_norm = [ones(ntr,1),(X - mean(X))./std(X)];
        
        Y = logit(chosen_dirdec(tr_in_slice,t));
        Y_norm = (Y - mean(Y))./std(Y);
        
        [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, Y_norm);
        beta(S,t,1:nX) = B_hat(2:end);
        cpd(S,t,1:nX) = CPD(2:end);
        ifsig(S,t,1:nX) = coeff_p(2:end)<thresh & model_p<thresh;
        pval(S,t,1:nX) = coeff_p(2:end);
    end
    
end

%% show regression results

ifsig(ifsig==0) = NaN;

figure;


for S = 1:nsubj
    
    % restrict trials
    subject = subject_names{S};
    ll = [];
    for T = 1:nterms
        
        [~,peak] = max(beta(S,:,T));
        disp([subject,': postprob ~ 1 + maxvalue'])
        disp(['    @peak=',num2str(1+t_mids(peak)),...
            ', B=',num2str(beta(S,peak,T)),...
            ', %CPD=',num2str(100*cpd(S,peak,T)),...
            ', p=',num2str(pval(S,peak,T))])
        
        subplot(2,2,S); hold on
        plot(t_mids,beta(S,:,T),'Color',[.5 .5 .5],'LineWidth',1);
        aa = plot(t_mids,ifsig(S,:,T).*beta(S,:,T),'Color',clrs_glm(T,:),'LineWidth',2);
        ll(T) = aa;
        
        if T==1
            xlim([-500 1500])
            ylim([-.05 .45])
            set(gca,'XTick',-500:500:1500,'XTickLabel',-.5:.5:1.5,...
                'YTick',0:.2:.4)
            plot([0 0 ],[-1 1],'k:')
            plot([-2000 2000],[0 0],'k:')
            
            xlabel('time from pics (s)')
            ylabel('\beta coefficient')
            
        end
        
        subplot(2,2,2+S); hold on
        plot(t_mids,100*cpd(S,:,T),'Color',[.5 .5 .5],'LineWidth',1);
        plot(t_mids,100*ifsig(S,:,T).*cpd(S,:,T),'Color',clrs_glm(T,:),'LineWidth',2);
        
        if T==1
            xlim([-500 1500])
            ylim([0 15])
            
             set(gca,'XTick',-500:500:1500,'XTickLabel',-.5:.5:1.5,...
                'YTick',0:5:15)
            
            plot([0 0 ],[0 100],'k:')
            
            xlabel('time from pics (s)')
            ylabel('%CPD')
        end
    end
    
    subplot(2,2,S)
%     legend(ll,glmvarnames,'box','off','Location','northeast')
end

%% check out gaze direction

if add_eyes
    
    [~,t0] = min(abs(t_eyes-0));
    
    figure;
    
    for S = 1:nsubj
        
        % restrict trials
        subject = subject_names{S};
        idx = find(strcmp(bhvdata.subject,subject_names{S}) & ...
            keep_tr);
        
        ntr = length(idx);
        
        rt = bhvdata.rt(idx);
        [rt_sort,row_sort] = sort(rt);
        medRT = median(rt);
        
        % if looking at chosen pic?
        gaze = bhvdata.location_target(idx,:);
        ll = bhvdata.lever(idx)==-1;
        gaze(ll,:) = -gaze(ll,:);
        for tr = 1:ntr
            bound = min(round(rt(tr)+2000+juice(S)),4000);
            gaze(tr,bound:end) = -2;
        end
        
        % plot heatmap
        subplot(4,2,S+[0,2,4]);
        imagesc(gaze(row_sort,:))
        hold on
        plot(t0+rt_sort,1:ntr,'k','LineWidth',2)
        
        % format
        plot([t0 t0],[1,ntr],'k','LineWidth',2)
        
        xlabel('time from pics (ms)')
        ylabel('trials')
        title(subject)
        
        set(gca,'XTick',1:500:4001,'XTickLabel',t_eyes(1:500:4001))
        
        xlim([1300,3500])
        
        
        
        % plot summary
        subplot(4,2,S+6);hold on
        missing = sum(gaze~=-2);
        
        plot(t_eyes,100*sum(gaze==1)./missing,'b','LineWidth',2);
        plot(t_eyes,100*sum(gaze==-1)./missing,'r','LineWidth',2);
        plot(t_eyes,100*sum(abs(gaze)==1)./missing,'m','LineWidth',1)
        
        % format
        plot([0 0],[0 100],'k','LineWidth',2)
        plot([-750 -750],[0 100],'k:','LineWidth',1)
        
        
        plot([medRT medRT],[0 100],'k:','LineWidth',2)
        xlabel('time from pics (ms)')
        ylabel('% trials')
        
        legend({'Ch','Unch','Ch+Unch'},'box','off','Location','northwest')
        
        xlim([-700,1500])
    end
    
    colormap([.8 .8 .8; 1 0 0 ; 1 1 1; 0 0 1])
    
end


end

