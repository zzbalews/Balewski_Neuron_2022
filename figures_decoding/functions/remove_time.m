function remove_time(bhvdata,valdata,keep_tr)
%%
% useful vars
subject_names = unique(valdata.subject);
nsubj = length(subject_names);

t_mids = valdata.t_mids;
t_subset = find(t_mids>=-500 & t_mids<1000);

figure;
for s = 1:nsubj
    
   % restrict to useful trials
   subject = subject_names{s};
   idx = strcmp(bhvdata.subject,subject) & keep_tr;
   ntr = sum(idx);
   
   %% plot avg ppd
   subplot(3,2,s); hold on
   
   shadedErrorBar(t_mids(t_subset),...
       mean(valdata.ch_ppd(idx,t_subset)),...
       std(valdata.ch_ppd(idx,t_subset))./sqrt(ntr),...
       'lineprops',{'Color',[1 0 0]});
   
    shadedErrorBar(t_mids(t_subset),...
       mean(valdata.unch_ppd(idx,t_subset)),...
       std(valdata.unch_ppd(idx,t_subset))./sqrt(ntr),...
       'lineprops',{'Color',[0 0 1]});
 
   YLIM = get(gca,'YLim');
   ylim(YLIM);
   xlim([-500 1000])
 plot([-2000 2000],[0 0],'k')
   plot([0 0],[-200 200],'k')
   xlabel('time from pics (ms)')
   ylabel({'OFC value decoding','%\Delta posterior probability'})
   title(subject)
   
   %% plot avg ppd residual, after linear time removed
   
   CH = valdata.ch_ppd(idx,t_subset);
   CHflat = reshape(CH,[],1);
   CH_m = mean(CHflat);
   CH_sd = std(CHflat);
   CHflat_norm = (CHflat - CH_m)./CH_sd;
   
   UNCH = valdata.unch_ppd(idx,t_subset);
   UNCHflat = reshape(UNCH,[],1);
   UNCH_m = mean(CHflat);
   UNCH_sd = std(UNCHflat);
   UNCHflat_norm = (UNCHflat - UNCH_m)./UNCH_sd;
   
   X = repmat(t_mids(t_subset),ntr,1);
   Xflat = reshape(X,[],1);
   temp = [Xflat];
   X_norm = [ones(size(Xflat)),(temp - mean(temp))./std(temp)];
   
   [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, CHflat_norm);
   CHflat_fit = CH_sd*(X_norm*B_hat) + CH_m;
   CHresid = reshape(CHflat - CHflat_fit,ntr,[]);
   
   [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, UNCHflat_norm);
   UNCHflat_fit = UNCH_sd*(X_norm*B_hat) + UNCH_m;
   UNCHresid = reshape(UNCHflat - UNCHflat_fit,ntr,[]);
   
   
   temp = reshape(CHflat_fit,ntr,[]);
   plot(t_mids(t_subset), mean(temp),'r--')
   temp = reshape(UNCHflat_fit,ntr,[]);
   plot(t_mids(t_subset), mean(temp),'b--')
   
   
   subplot(3,2,s+2); hold on
   shadedErrorBar(t_mids(t_subset), ...
       mean(CHresid),...
       std(CHresid)./sqrt(ntr),...
       'lineprops',{'Color',[1 0 0]})
   shadedErrorBar(t_mids(t_subset), ...
       mean(UNCHresid),...
       std(UNCHresid)./sqrt(ntr),...
       'lineprops',{'Color',[0 0 1]})
   
   YLIM = get(gca,'YLim');
   ylim(YLIM);
   plot([0 0],[-200 200],'k')
   plot([-2000 2000],[0 0],'k')
    
   xlim([-500 1000])
   xlabel('time from pics (ms)')
    ylabel('residual %\DeltaPP, after linear time')
    
    %% plot avg ppd, after remove mean
    subplot(3,2,s+4); hold on
    CH_m = repmat(mean(valdata.ch_ppd(idx,t_subset)),ntr,1);
    UNCH_m = repmat(mean(valdata.unch_ppd(idx,t_subset)),ntr,1);
    
    CH_0 = valdata.ch_ppd(idx,t_subset) - CH_m;
    UNCH_0 = valdata.unch_ppd(idx,t_subset) - UNCH_m;
    
    shadedErrorBar(t_mids(t_subset),...
        mean(CH_0),...
        std(CH_0)/sqrt(ntr),...
        'lineprops',{'Color',[1 0 0]})
    
    shadedErrorBar(t_mids(t_subset),...
        mean(UNCH_0),...
        std(UNCH_0)/sqrt(ntr),...
        'lineprops',{'Color',[0 0 1]})
   
    
     YLIM = get(gca,'YLim');
   ylim(YLIM);
   plot([0 0],[-200 200],'k')
   plot([-2000 2000],[0 0],'k')
    
   xlim([-500 1000])
   
    xlabel('time from pics (ms)')
    ylabel('%\DeltaPP - mean')
   
end



end

