function do_free_apply(minidata,pl2dir,spksmooth,spkslice,outdir,decoder,niters)

%% data details

temp = strsplit(spksmooth{1},'_');
session = strjoin(temp(1:3),' ');
if length(spksmooth)==1
    region = strrep(temp{end},'.mat','');
else
    region = strrep(temp{end}(2:end),'.mat','');
end


% outprefix_old = strjoin({[outdir,'/old'],[strrep(session,' ','_'),'_',region,'_value_free_w400']},'/');
outprefix = strjoin({outdir,[strrep(session,' ','_'),'_',region,'_value_free_w400']},'/');

% set up decoder details
if strcmp(decoder, 'value')
    levels = 1:4;
    nlevels = 4;
    
    truelabel = minidata.valbin;
    
else
    disp('need to set up this decoder! STOP')
    return
end
% 
%% select trials

% use forced to train decoder
[~,temp] = max(minidata.valbin,[],2);
highvalloc = temp*2-3;

forced_tr = minidata.trialtype==1 & minidata.lever==highvalloc; %minidata.lever~=0;% 

% make balanced sets for training
[trialsLOO, setsize, N] = pull_balanced_value(minidata, forced_tr, niters);

% apply to free 
free_tr = minidata.trialtype==2 & minidata.lever~=0;
ntr = sum(free_tr);

%% load brain data
varnames = {'SPKfrnorm_units',...
    'LFP01delta','LFP02theta','LFP03alpha',...
    'LFP04beta','LFP05gamma','LFP06highgamma'};

% varnames = {'SPKfrnorm_units'};

brain_slice = cell(1,length(varnames));
brain_smooth = cell(1,length(varnames));

unit_names = {};
for i = 1:length(spksmooth)
    
    % load slice (forced; train)
    load(strjoin({pl2dir,spksmooth{i}},'/'),...
        'SPKfrnorm_names','t_mids');
    unit_names = cat(1,unit_names,SPKfrnorm_names);
    
    for v = 1:length(varnames)
        slice = load(strjoin({pl2dir,spkslice{i}},'/'), varnames{v});
        brain_slice{v} = cat(3, brain_slice{v}, slice.(varnames{v}));
        
        smooth = load(strjoin({pl2dir,spksmooth{i}},'/'), varnames{v});
        brain_smooth{v} = cat(3, brain_smooth{v}, smooth.(varnames{v}));
    end
end

nmids = length(t_mids);

%% cylce through iters, average 
PCthresh = 95;

free_val = truelabel(free_tr,:);
free_postprob = nan(ntr,nmids,nlevels,niters);

for i = 1:niters
    tic
    % get forced trials in this iter
    input_tr = reshape(trialsLOO{i,1},[],1);
    input_label = max(truelabel(input_tr,:),[],2);
    
    % get top pcs for slice
    input_pcs = [];
    apply_pcs = [];
    for freq = 1:length(varnames)
        snippet_train = permute(brain_slice{freq}(input_tr,:,:),[1,3,2]);
        [coeff, score, ~, ~, expvar, mu] = pca(snippet_train);
        PC95 = find(cumsum(expvar)>PCthresh,1);
        
        input_pcs = cat(2, input_pcs, score(:,1:PC95));

    % use same pc projections for smooth
    snippet_apply = permute(brain_smooth{freq}(free_tr,:,:),[1,3,2]);
    
    term1 = snippet_apply - repmat(mu,ntr,1,nmids);
    term2 = repmat(coeff(:,1:PC95),1,1,nmids);
    
% % % %     apply_pcs = nan(ntr,PC95,nmids);
% % % %     for m = 1:nmids
% % % %         apply_pcs(:,:,m) = term1(:,:,m) * term2(:,:,m);
% % % %     end
    apply_pcs = cat(2, apply_pcs, pagemtimes(term1,term2)); % ntr x pcs x nmids
    end
    
    % train decoder
    mdl_iter = fitcdiscr(input_pcs,input_label); % use default LDA
    
    % apply decoder to all times in smooth
% t_idx = find(t_mids>-1600 & t_mids<2000);
% t_idx = find(t_mids>-1100 & t_mids<1100);

t_idx = find(t_mids>-1600 & t_mids<1100);

    for m = t_idx % 1:nmids
        [apply_guess, apply_postprob, ~] = predict(mdl_iter,apply_pcs(:,:,m));
        free_postprob(:,m,:,i) = apply_postprob;
    end
    toc
    
end

% average across iters
free_postprob_avg = mean(free_postprob,4);

%% get NA baselines

% find NA values for each trial
na_vals = nan(ntr,3);
for tr = 1:ntr
    
   temp = levels(~ismember(levels,free_val(tr,:))); 
   na_vals(tr,1:length(temp)) = temp;
   
end

% calculate NA baseline
baseline_postprob = nan(1,nmids,nlevels);
for L = 1:nlevels
    % find trials where value is missing
    na_trs = sum(na_vals==levels(L),2)>0;
        baseline_postprob(:,:,L) = mean(free_postprob_avg(na_trs,:,L));
end

% figure; hold on
% 
% for i = 1:4
%     plot(t_mids,baseline_postprob(:,:,i))
% end
% plot([-2000 2000],[.25 .25],'k:')
% plot([0,0],[0,.5],'k')
% ylim([0 .5])
% legend({'v=1','v=2','v=3','v=4'},'box','off')
% xlabel('time from pics (ms)')
% ylabel('NA baseline')
% temp = strsplit(outprefix,'/');
% 
% title(strrep(temp{end},'_','\_'))

%% get postprob %delta from baseline

free_postprob_delta = 100*(free_postprob_avg - baseline_postprob) ./ baseline_postprob;

% % %% find chosen/unchosen post prob delta
% % ch_ppd = nan(ntr,nmids);
% % unch_ppd = nan(ntr,nmids);
% % 
% % for tr = 1:ntr
% %     ch = max(free_val(tr,:));
% %     unch = min(free_val(tr,:));
% %     ch_ppd(tr,:) = free_postprob_delta(tr,:,ch);
% %     unch_ppd(tr,:) = free_postprob_delta(tr,:,unch);
% % 
% % end
% % 
% % save([outprefix,'.mat'],'free_postprob_delta','free_tr','ch_ppd','unch_ppd','t_mids');

% %% load saved decoding
% load([outprefix_old,'.mat'])

%% get states
thresh = 200;
free_postprob_delta_thresh = free_postprob_delta>=200;

[ntr,nmids,nlevels] = size(free_postprob_delta_thresh);
free_postprob_states = nan(ntr,nmids,nlevels);

for i = 1:ntr
    for j = 1:nlevels

        temp = free_postprob_delta_thresh(i,:,j);
        free_postprob_states(i,:,j) = give_consec_seg(temp,4);

    end
end


%% relabel states by ch/unch/na1/na2

free_val = minidata.valbin(free_tr,:);
lever = minidata.lever(free_tr);
lever(lever==0) = NaN;
lever = lever/2 + 1.5;

levels = 1:4;

ch_ppd = nan(ntr,nmids);
unch_ppd = nan(ntr,nmids);
na_ppd = nan(ntr,nmids,3);

ch_state = nan(ntr,nmids);
unch_state = nan(ntr,nmids);
na_state = nan(ntr,nmids,3);

for i = 1:ntr
   
    ch = free_val(i,lever(i));
    unch = free_val(i,3-lever(i));
    
    na = levels(~ismember(levels,[ch,unch]));
    
    ch_ppd(i,:) = free_postprob_delta(i,:,ch);
    unch_ppd(i,:) = free_postprob_delta(i,:,unch);
    na_ppd(i,:,1:length(na)) = free_postprob_delta(i,:,na);

    ch_state(i,:) = free_postprob_states(i,:,ch);
    unch_state(i,:) = free_postprob_states(i,:,unch);
    na_state(i,:,1:length(na)) = free_postprob_states(i,:,na);
    
end


%% save new files
 save([outprefix,'.mat'],'free_postprob_delta','free_tr','t_mids',...
     'ch_ppd','unch_ppd','na_ppd','ch_state','unch_state','na_state');

%% show avg ppd_delta
figure; hold on

plot([-2000 2000],[.25 .25],'k:')
plot([0 0],[-50 150],'k')
plot(t_mids, mean(ch_ppd),'LineWidth',2,'Color',[.9 0 0])
plot(t_mids, mean(unch_ppd),'LineWidth',2,'Color',[0 0 .9])

xlim([-800 1500])
ylim([-20 120])

xlabel('time from pics (ms)')
ylabel('%\Delta posterior probability value decoder')

print([outprefix,'.png'],'-dpng')
end
