function checkout_decoder_latency(bhvdata,valdataOFC,dirdata,keep_tr)

% useful vars


subject_names = unique(bhvdata.subject);
nsubj = length(subject_names);

t_mids = dirdata.t_mids;
[~,t0] = min(abs(t_mids-0));
nmids = length(t_mids);

ntr = sum(keep_tr);

%% get dirdata 1st state latency (option 1 = ignore ch/unch, option 2 = only ch)

CNstates = dirdata.postprob_ch_disc==1;% | dirdata.postprob_ch_disc==-1;
CNstates = CNstates(keep_tr,:);

CNstates_onset = nan(ntr,1);
for tr = 1:ntr
   
    idx = find(CNstates(tr,t0:end),1,'first');
    if ~isempty(idx)
    CNstates_onset(tr) = t_mids(t0+idx-1);
    end
    
end


%% get valdataOFC 1st state & latency
OFCstates = valdataOFC.ch_state;
% OFCstates = valdataOFC.ch_state - valdataOFC.unch_state;
OFCstates = OFCstates(keep_tr,:);

OFCstates_onset = nan(ntr,2);
[~,t100] = min(abs(t_mids-50));
for tr = 1:ntr

    idx = find(OFCstates(tr,t100:end)~=0,1,'first');
    if ~isempty(idx)
       OFCstates_onset(tr,1) = t_mids(t100+idx-1);
       OFCstates_onset(tr,2) = OFCstates(tr,t100+idx-1);
    end
    
end


%%

figure; 
for s = 1:2
    for v = 1:3
    subplot(3,2,(v-1)*2+s); hold on
    idx = strcmp(bhvdata.subject(keep_tr),subject_names{s}) & (v+1)==max(bhvdata.valbin_expval(keep_tr,:),[],2);
    histogram(CNstates_onset(idx),'BinWidth',20,'Normalization','probability')

    histogram(OFCstates_onset(idx,1),'BinWidth',20,'Normalization','probability')
    
    
    xlim([0 1000])
    xlabel('time from pics (ms)')
    ylabel('norm. trial count')
    end
end

legend({'CN','OFC'})

%% compare 1st CN state latency with strength of OFC value decoder


track_B = nan(2,nmids,s);
track_cpd = nan(2,nmids,s);
track_coeff_p = nan(2,nmids,s);
track_model_p = nan(1,nmids,s);
for s = 1:2
X = [CNstates_onset, max(bhvdata.valbin_expval(keep_tr,:),[],2)];
oktrs = ~isnan(X(:,1)) & strcmp(bhvdata.subject(keep_tr),subject_names{s});
X = X(oktrs,:);
% vif = max(diag(inv(corrcoef(X))));

X_norm = [ones(sum(oktrs),1), (X - mean(X))./std(X)];

Y = valdataOFC.unch_ppd(keep_tr,:);
Y = Y(oktrs,:);
Y_norm = (Y - mean(Y))./std(Y);

for t = 1:nmids
    [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, Y_norm(:,t));
    
    track_B(:,t,s) = B_hat(2:end);
    track_cpd(:,t,s) = CPD(2:end);
    track_coeff_p(:,t,s) = coeff_p(2:end);
    track_model_p(1,t,s) = model_p;
end

end
figure; 
thresh = 0.0001;
ifsig = double(track_model_p < thresh & track_coeff_p < thresh);
ifsig(ifsig==0) = NaN;

for s = 1:2
subplot(2,2,s); hold on
plot(t_mids, track_B(:,:,s),'Color',[.5 .5 .5],'LineWidth',1)
p = plot(t_mids, track_B(:,:,s) .* ifsig(:,:,s),'LineWidth',2);
plot([-2000 2000],[0 0],'k:')
plot([0 0],[-.3 .3],'k:')
ylim([-.25 .1])
ylabel('\beta coefficient')
legend([p(1),p(2)],{'CN state onset','maxval'},'box','off','Location','southwest')

subplot(2,2,2+s); hold on
plot(t_mids, 100*track_cpd(:,:,s),'Color',[.5 .5 .5],'LineWidth',1)
plot(t_mids, 100*track_cpd(:,:,s) .* ifsig(:,:,s),'LineWidth',2)
plot([0 0],[0 5],'k:')
ylim([0 5])
ylabel('% CPD')

end
for i =1:4
    subplot(2,2,i)
    xlim([-750 1000])
    xlabel('time from pics (ms)')
end

%% compare which OFC state is first following CN direction onset
OFCstates = valdataOFC.ch_state - valdataOFC.unch_state;
OFCstates = OFCstates(keep_tr,:);

OFCstates_onset = nan(ntr,2);
[~,t100] = min(abs(t_mids-50));
track = nan(ntr,1);
for tr = 1:ntr

    if isnan(CNstates_onset(tr))
        continue
    end
    [~,t_start] = min(abs(t_mids - CNstates_onset(tr)));
    track(tr) = t_start;
    idx = find(OFCstates(tr,t_start:end)~=0,1,'first');
    if ~isempty(idx)
       OFCstates_onset(tr,1) = t_mids(t_start+idx-1);
       OFCstates_onset(tr,2) = OFCstates(tr,t_start+idx-1);
    end
    
end

figure; hold on
histogram(CNstates_onset,'BinWidth',20)
histogram(OFCstates_onset(:,1),'BinWidth',20)

figure; 
temp = OFCstates_onset(:,1)-CNstates_onset;
histogram(temp);

figure; 
histogram(OFCstates_onset(:,2))



%% split trials by if chosen CN dir = 1st state --> looks the same as split by 1st saccade

niters = 1000;
tr_sampled = sample_matched_trials(bhvdata,keep_tr,niters,0);
    

temp = dirdata.postprob_ch_disc==1;
CN1stch = nan(size(temp,1),1);
for tr = 1:size(temp,1)
   idx = find(temp(tr,t0:end)==1,1,'first');
   if ~isempty(idx)
       CN1stch(tr) = 1;
   end
end

clrs = [1 .6 0;.2 .2 .2];

pp = {valdataOFC.ch_ppd,valdataOFC.unch_ppd};
XLIM = [-500 1000];

YLIM = [-50 200];
    YTICK = 0:100:200;
    yy = {'chosen','unchosen'};
    decoder = 'value';
        metric = '%\Delta posterior probability';

figure;
for S = 1:nsubj
    subject = subject_names{S};
    idx = keep_tr & strcmp(bhvdata.subject,subject);
    
    
    for r = 1:2
        % plot avg post prob
        subplot(2,2,S+2*(r-1)); hold on
        plot(XLIM,[0 0],'k:')
        plot([0 0],YLIM,'k:')
        
        temp = nan(niters,nmids);
        for j = 1:nmids
            HH = reshape(tr_sampled.(subject).keep_match,1,[]);
            KK = pp{r}(HH,j)';
            LL = reshape(KK,[],niters);
            temp(:,j) = mean(LL);
        end

        Y = quantile(temp,[.001,.5,.999],1);
        shadedErrorBar(t_mids, Y(2,:),...
            [Y(1,:)-Y(2,:);Y(2,:)-Y(3,:)],'lineprops',{'Color',clrs(1,:)});
        
        a = plot([NaN NaN],[.5 .5],'Color',clrs(1,:));
        b = plot(t_mids, mean(pp{r}(tr_sampled.(subject).keep_mis,:)),'Color',clrs(2,:),'LineWidth',1);
  
        % format
        text(-250,mean(YLIM),['n=',num2str(length(tr_sampled.(subject).keep_mis))],'HorizontalAlignment','center')
        xlim(XLIM)
        ylim(YLIM)
        
        set(gca,'XTick',-500:500:1500,'XTickLabel',-.5:.5:1.5,...
            'YTick',YTICK)

        xlabel('Time from pictures (s)')
        ylabel([yy{r},' ',decoder,' decoder: ',metric])
    end
end




end