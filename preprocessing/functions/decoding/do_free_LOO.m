function [outfile_prefix] = do_free_LOO(minidata,pl2dir,spksmooth,outdir,decoder,niters)

%% data details
temp = strsplit(spksmooth{1},'_');
session = strjoin(temp(1:3),' ');
region = strrep(temp{end},'.mat','');

outfile_prefix = strjoin({outdir,[strrep(session,' ','_'),'_',region,'_',decoder,'_freeLOO']},'/');
%% set up decoder details
if strcmp(decoder,'direction')
    
    levels = [-1,1];
    nlevels = 2;
    
    truelabel = minidata.lever;
elseif strcmp(decoder,'sacc1') 
    
    levels = [-1,1];
    nlevels = 2;
    
    truelabel = minidata.saccloc(:,1);

elseif strcmp(decoder,'sacc2')
    
    levels = [-1,1];
    nlevels = 2;
    
    truelabel = minidata.saccloc(:,2);


else
    disp('need to set up this decoder! STOP')
    return
end

%% select trials
keep_tr = minidata.trialtype==2 & minidata.lever~=0;
% keep_tr = minidata.trialtype==2 & minidata.lever~=0 & ~isnan(minidata.sacc_to_lever(:,2));

% make balanced sets
[trials_balanced, setsize, N, leftover] = pull_balanced_direction(minidata, truelabel, keep_tr, niters);
setsize,N
trialsLOO = trials_balanced(:,1); % balanced sets to use for LOO
trialsREM = trials_balanced(:,2); % all other free trials

ntr = (setsize + leftover) * N;
%% load brain data

% % %
% % % braindata = [];
% % %
% % % for f = 1:length(spksmooth)
% % %     load(strjoin({pl2dir,spksmooth{f}},'/'),...
% % %         't_mids','SPKfrnorm_units');
% % %
% % %     braindata = cat(3,braindata,SPKfrnorm_units);
% % % end



% 
% varnames = {'SPKfrnorm_units',...
%     'LFP01delta','LFP02theta','LFP03alpha',...
%     'LFP04beta','LFP05gamma','LFP06highgamma'};

varnames = {'SPKfrnorm_units'};


braindata = cell(1,length(varnames));
for f = 1:length(spksmooth)
    load(strjoin({pl2dir,spksmooth{f}},'/'),'t_mids');
    
    
    for v = 1:length(varnames)
        data = load(strjoin({pl2dir,spksmooth{f}},'/'),varnames{v});
        braindata{v} = cat(3,braindata{v},data.(varnames{v}));
        
        clear(varnames{v})
    end
end

nmids = length(t_mids);

disp('data loaded!')
%% set up output structs
LDA_output = struct();
LDA_output.trialgroup = repelem({'LOO','remaining'},N*[setsize, leftover]);

% set up arrays to collect parfor output
collect_trialnumber =  nan(ntr, nmids, niters);
collect_truelabel = nan(ntr, nmids, niters);
collect_guess = nan(ntr, nmids, niters);
collect_guess_shuffle = nan(ntr, nmids, niters);
collect_postprob = nan(ntr, nmids, niters, nlevels);
collect_postprob_shuffle = nan(ntr, nmids, niters, nlevels);

%% start decoding!

t_window = find(t_mids>=-1100 & t_mids<=1600);
% t_window = find(t_mids>=-1600 & t_mids<=1100);

parfor m = t_window %1:nmids% 249:749%48:nmids
    if mod(m,5)==0
        disp(['... t = ',num2str(t_mids(m))])
    end
    
    % % %     % use all unit firing rates
    % % %     braindata_pc = permute(braindata(:,m,:),[1,3,2]);
    
    
    % do PCA on this slice of brain data & trials; keep top 95% var PCs for each freq band & spk
    braindata_pc = [];
    for freq = 1:length(braindata)
        
        [~,score,~,~,expvar] = pca(permute(braindata{freq}(:,m,:),[1,3,2]));
        PC95 = find(cumsum(expvar)>95,1);
        braindata_pc = cat(2,braindata_pc,score(:,1:PC95));
        
    end
    
    % set up arrays for slice outputs
    slice_trialnumber = nan(ntr, 1, niters);
    slice_truelabel = nan(ntr, 1, niters);
    slice_guess = nan(ntr, 1, niters);
    slice_guess_shuffle = nan(ntr, 1, niters);
    slice_postprob = nan(ntr, 1, niters, nlevels);
    slice_postprob_shuffle = nan(ntr, 1, niters, nlevels);
    
    for i = 1:niters % for each shuffled trial iteration)
        for s = 1:N % for each LOO set
            
            % split train/test trials
            tr_train = reshape(trialsLOO{i}(:,[1:s-1, s+1:end]),[],1);
            tr_test = trialsLOO{i}(:,s);
            tr_rem = trialsREM{i}(:,s);
            tr_rem = tr_rem(~isnan(tr_rem));
            
            % get brain data for trials
            input_train = braindata_pc(tr_train,:);
            input_test = braindata_pc(tr_test,:);
            input_rem = braindata_pc(tr_rem,:);
            
            % get correct labels for trials, plus shuffle
            label_train = truelabel(tr_train);
            label_train_shuffle = label_train(randperm(length(label_train)));
            label_true = truelabel(tr_test);
            label_rem = truelabel(tr_rem);
            
            % do LDA
            try
                [guess_test, ~, postprob_test] = classify(input_test, input_train, ...
                    label_train,'diaglinear','empirical'); % true LOO
            catch
                guess_test = nan(size(label_true));
                postprob_test = nan(length(label_true),nlevels);
            end
            
            try
                [guess_test_shuffle, ~, postprob_test_shuffle] = classify(input_test, input_train, ...
                    label_train_shuffle,'diaglinear','empirical'); % shuffle training labels for LOO
            catch
                guess_test_shuffle = nan(size(label_true));
                postprob_test_shuffle = nan(length(label_true),nlevels);
            end
            
            try
                [guess_rem, ~, postprob_rem] = classify(input_rem, input_train, ...
                    label_train,'diaglinear','empirical'); % true LOO %%% diaglinear
            catch
                guess_rem = nan(size(label_rem));
                postprob_rem = nan(length(label_rem),nlevels);
            end
            
            % update slice arrays
            
            idx_test = (s-1)*setsize + (1:setsize)';
            idx_rem = (s-1)*leftover + (1:length(tr_rem))' + setsize*N;
            
            slice_trialnumber([idx_test;idx_rem],1,i) = [tr_test;tr_rem];
            
            slice_truelabel([idx_test;idx_rem],1,i) = [label_true;label_rem];
            
            slice_guess([idx_test;idx_rem],1,i) = [guess_test;guess_rem];
            slice_guess_shuffle(idx_test,1,i) = guess_test_shuffle;
            
            slice_postprob([idx_test;idx_rem],1,i,:) = [postprob_test;postprob_rem];
            slice_postprob_shuffle(idx_test,1,i,:) = postprob_test_shuffle;
        end
    end
    
    % save slice output
    collect_trialnumber(:,m,:) = slice_trialnumber;
    
    collect_truelabel(:,m,:) = slice_truelabel;
    
    collect_guess(:,m,:) = slice_guess;
    collect_guess_shuffle(:,m,:) = slice_guess_shuffle;
    
    collect_postprob(:,m,:,:) = slice_postprob;
    collect_postprob_shuffle(:,m,:,:) = slice_postprob_shuffle;
    
end

% update output from all time slices
LDA_output.trialnumber = collect_trialnumber;
LDA_output.truelabel = collect_truelabel;
LDA_output.guess = collect_guess;
LDA_output.guess_shuffle = collect_guess_shuffle;
LDA_output.postprob = collect_postprob;
LDA_output.postprob_shuffle = collect_postprob_shuffle;



%% get states from post prob


% re-organize postprobs from all iters by trial #
alltrials_postprob = nan(size(braindata{1},1), nmids, 1, 2, niters); % trials x time x trialperiod x levels x niters

idx = find(t_mids==-1);

for i = 1:niters
    
    tr_current = LDA_output.trialnumber(:,idx,i);
    keep = ~isnan(tr_current);
    
    alltrials_postprob(tr_current(keep),:,:,:,i) = LDA_output.postprob(keep,:,i,:);
    
end

% avg across all iters
alltrials_postprob_avg = nanmean(alltrials_postprob,5);

% get states
thresh = [.75,.8,.85,.9,.95];
alltrials_states = get_states(alltrials_postprob_avg, 2, 4, thresh);

%% save outputs

save([outfile_prefix,'.mat'],'LDA_output','alltrials_postprob_avg','alltrials_states','t_mids','thresh','-v7.3')


%% show LOO accuracy

tr_LOO = strcmp(LDA_output.trialgroup,'LOO');

figure; hold on

plot([0 0],[0 1],'k:')
m_shuf = nanmean(nanmean(LDA_output.guess_shuffle(tr_LOO,:,:)==LDA_output.truelabel(tr_LOO,:,:),3),1);
plot(t_mids,m_shuf,'Color',[.5 .5 .5])

m = nanmean(nanmean(LDA_output.guess(tr_LOO,:,:)==LDA_output.truelabel(tr_LOO,:,:),3),1);
plot(t_mids,m,'Color',[0 .7 0],'LineWidth',2)


xlabel('time from pics (ms)')
ylabel(['LOO accuracy (',decoder,', free)'])

% xlim([-1500 1000])
xlim([-1000 1500])
ylim([.4 1])

print([outfile_prefix,'.png'],'-dpng')
end

