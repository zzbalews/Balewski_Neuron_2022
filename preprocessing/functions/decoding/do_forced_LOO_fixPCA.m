function do_forced_LOO_fixPCA(minidata,pl2dir,spksmooth,outdir,decoder,niters,trialperiod)

%% data details
temp = strsplit(spksmooth{1},'_');
session = strjoin(temp(1:3),' ');
if length(spksmooth)==1
    region = strrep(temp{end},'.mat','');
else
    region = strrep(temp{end}(2:end),'.mat','');
end

out_prefix = strjoin({outdir,[strrep(session,' ','_'),'_',region,'_value_forcedLOO']},'/');
%% set up decoder details
if strcmp(decoder,'value')
    
    levels = 1:4;
    nlevels = 4;
    
    truelabel = max(minidata.valbin,[],2);
    
else
    disp('need to set up this decoder! STOP')
    return
end

%% select trials
% limit to forced
[~,temp] = max(minidata.valbin,[],2);
highvalloc = temp*2-3;

keep_tr = minidata.trialtype==1 & minidata.lever==highvalloc; %minidata.lever~=0;% minidata.lever==highvalloc; %

% make balanced sets
[trialsLOO, setsize, N] = pull_balanced_value(minidata, keep_tr, niters);
ntr = length(reshape(trialsLOO{1,1},[],1));

%% load brain data
% varnames = {'SPKfrnorm_units',...
%     'LFP01delta','LFP02theta','LFP03alpha',...
%     'LFP04beta','LFP05gamma','LFP06highgamma'};
% 
varnames = {'SPKfrnorm_units'};

braindata = cell(1,length(varnames));

unit_names = {};
for i = 1:length(spksmooth)
    load(strjoin({pl2dir,spksmooth{i}},'/'),'t_mids','SPKfrnorm_names');
    unit_names = cat(1,unit_names,SPKfrnorm_names);
    for v = 1:length(varnames)
        data = load(strjoin({pl2dir,spksmooth{i}},'/'),varnames{v});
        braindata{v} = cat(3, braindata{v}, data.(varnames{v}));
        
        clear('data');
    end
end
nmids = length(t_mids);

disp('data loaded!')

% %% restrict units
% unit_frs = nan(size(unit_names));
% unit_quals = nan(size(unit_names));
% for u = 1:length(unit_frs)
%    
%     idx = find(strcmp(unit_names{u},unit_info.name));
%     if ~isempty(idx)
%     unit_frs(u) = unit_info.fr(idx);
%     unit_quals(u) = unit_info.qual(idx);
% end
% end
% 
% % keep_units = unit_quals>2;
% % braindata{1} = braindata{1}(:,:,keep_units);


%% cycle through sets (LOO), and iters

decoder_guess = nan(ntr, nmids, niters);
decoder_shuf = nan(ntr,nmids,niters);
decoder_tru = nan(ntr,nmids, niters);

if strcmp(trialperiod,'pics')
t_idx = find(t_mids>=-800 & t_mids<=1550);
elseif strcmp(trialperiod,'choice')
    t_idx = find(t_mids>=-1550 & t_mids<=800);
end

PCthresh = 95;

parfor i = 1:niters
    tic
    
    iter_guess = nan(ntr, nmids);
    iter_shuf = nan(ntr,nmids);
    iter_tru = nan(ntr,nmids);
    for m = t_idx%1:nmids % time points
        
        
        for S = 1:N % cylce through LOO sets
            
             % trials in this set
        tr_set = reshape(trialsLOO{i,1},[],1);
        
            % split into test/train trials
            tr_test = (S-1)*setsize + (1:setsize)';
            temp = 1:length(tr_set);
            tr_train = find(~ismember(temp,tr_test))';
            
            ntr_test = length(tr_test);
            
       
        
        % load brain data; do PCA (train only), keep top 95% var PCs
        input_train = [];
        input_test = [];
        for freq = 1:length(braindata)
            
            % PCA on training set
            snippet_train = permute(braindata{freq}(tr_set(tr_train),m,:),[1,3,2]);
            [coeff, score, ~, ~, expvar, mu] = pca(snippet_train);
            PC95 = find(cumsum(expvar)>PCthresh,1);
            
            % update 
            term1 = snippet_train - repmat(mu, size(snippet_train,1),1);
            term2 = coeff(:,1:PC95);
            
            input_train = cat(2, input_train, term1*term2);
        
            % same PC projection on testing set
            snippet_test = permute(braindata{freq}(tr_set(tr_test),m,:),[1,3,2]);
            term1 = snippet_test - repmat(mu, size(snippet_test,1),1);
            term2 = coeff(:,1:PC95);
            
            input_test = cat(2, input_test, term1*term2);
            
        end
        
        
            % get correct labels
            label_test_tru = truelabel(tr_set(tr_test));
            label_train_tru = truelabel(tr_set(tr_train));
            label_train_shuf = label_train_tru(randperm(length(tr_train)));
           
            % do LDA
%            try
            [guess_test, ~, postprob_test] = classify(input_test, input_train, label_train_tru,...
                'linear','empirical'); % true LOO
%            catch
%                guess_test = nan(ntr_test,1);
%                postprob_test= nan(ntr_test,4);
%            end
               
%            try
            [guess_test_shuf, ~, postprob_test_shuf] = classify(input_test, input_train, label_train_shuf,...
                'linear','empirical'); % shuffle LOO
%              catch
%                guess_test_shuf = nan(ntr_test,1);
%                postprob_test_shuf= nan(ntr_test,4);
%            end
%            
            % save result
            iter_guess(tr_test,m) = guess_test;
            iter_shuf(tr_test,m) = guess_test_shuf;
            iter_tru(tr_test,m) = label_test_tru;
            
        end
        
    end
    decoder_guess(:,:,i) = iter_guess;
    decoder_shuf(:,:,i) = iter_shuf;
    decoder_tru(:,:,i) = iter_tru;
    
    toc
    
end

%% show accuracy
figure;
subplot(2,1,1); hold on
m = mean(nanmean(decoder_tru==decoder_guess,3),1);
plot(t_mids,m,'Color',[0 .7 0]);
m_shuf = mean(nanmean(decoder_tru==decoder_shuf,3),1);
plot(t_mids,m_shuf,'Color',[.5 .5 .5]);

plot([0 0],[0 1],'k')
ylim([.2 .6])

if strcmp(trialperiod,'pics')
xlim([-750,1500])
xlabel('time from pics (ms)')
elseif strcmp(trialperiod,'choice')
    xlim([-1500,750])
xlabel('time from choice (ms)')
end
ylabel('value decoder accuracy')

title({session,[region,' (PC exp var ',num2str(PCthresh),')']})

%% show conf mat at peak accuracy
[~,peak] = max(m);

if strcmp(trialperiod,'pics')
peak = min(peak, 48);
elseif strcmp(trialperiod,'choice')
    peak = min(peak,39);
end
% peak = 45;
peak_tru = reshape(decoder_tru(:,peak,:),[],1);
peak_guess = reshape(decoder_guess(:,peak,:),[],1);
peak_shuf = reshape(decoder_shuf(:,peak,:),[],1);

C = confusionmat(peak_tru,peak_guess,'order',levels);
C = C';

total = repmat(sum(C,1),nlevels,1);
conf_mat = C./total;

subplot(2,1,2);
imagesc(conf_mat)
daspect([1 1 1])
colorbar
colormap(hot)
caxis([.25 .8])
set(gca,'XTick',1:4,'XTickLabel',levels,...
    'YTick',1:4,'YTickLabel',levels)
xlabel('true')
ylabel('decoder')

title(['peak @ t=',num2str(t_mids(peak))])

%% format & save
set(gcf,'Position',[100 100 250 450])
print([out_prefix,'.png'],'-dpng')


save([out_prefix,'.mat'],'t_mids','m','conf_mat','-v7.3')


end

