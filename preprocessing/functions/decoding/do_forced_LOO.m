function do_forced_LOO(minidata,pl2dir,spksmooth,units,outdir,decoder,niters)

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
braindata = [];


varnames = {'SPKfrnorm_units',...
    'LFP01delta','LFP02theta','LFP03alpha',...
    'LFP04beta','LFP05gamma','LFP06highgamma'};

% varnames = {'SPKfrnorm_units'};

braindata = cell(1,length(varnames));


for i = 1:length(spksmooth)
    load(strjoin({pl2dir,spksmooth{i}},'/'),'t_mids');
    
    for v = 1:length(varnames)
        data = load(strjoin({pl2dir,spksmooth{i}},'/'),varnames{v});
        braindata{v} = cat(3, braindata{v}, data.(varnames{v}));
        
        clear(varnames{v});
    end
end
nmids = length(t_mids);

disp('data loaded!')
% %% restrict to units >1Hz
% keep_units = units.fr>=1;
% braindata = braindata(:,:,keep_units);


%% cycle through sets (LOO), and iters

decoder_guess = nan(ntr, nmids, niters);
decoder_shuf = nan(ntr,nmids,niters);
decoder_tru = nan(ntr,nmids, niters);

t_idx = find(t_mids>=-800 & t_mids<=1550);

PCthresh = 95;

parfor i = 1:niters
    tic
    
    iter_guess = nan(ntr, nmids);
    iter_shuf = nan(ntr,nmids);
    iter_tru = nan(ntr,nmids);
    for m = t_idx%1:nmids % time points
        
        % trials in this set
        tr_set = reshape(trialsLOO{i,1},[],1);
        
        % do PCA, keep top 95% var PCs
        braindata_pc = [];
        for freq = 1:length(braindata)
            
            snippet = permute(braindata{freq}(tr_set,m,:),[1,3,2]);
            [~,score,~,~,expvar] = pca(snippet);
            PC95 = find(cumsum(expvar)>PCthresh,1);
            braindata_pc = cat(2,braindata_pc,score(:,1:PC95));
        
        end
        
        
        for S = 1:N % cylce through LOO sets
            % split into test/train trials
            tr_test = (S-1)*setsize + (1:setsize)';
            temp = 1:length(tr_set);
            tr_train = find(~ismember(temp,tr_test))';
            
            % split brain data
            input_test = braindata_pc(tr_test,:);
            input_train = braindata_pc(tr_train,:);
            
            % get correct labels
            label_test_tru = truelabel(tr_set(tr_test));
            label_train_tru = truelabel(tr_set(tr_train));
            label_train_shuf = label_train_tru(randperm(length(tr_train)));
            
            % do LDA
            [guess_test, ~, postprob_test] = classify(input_test, input_train, label_train_tru,...
                'linear','empirical'); % true LOO
            
            [guess_test_shuf, ~, postprob_test] = classify(input_test, input_train, label_train_shuf,...
                'linear','empirical'); % shuffle LOO
            
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
m = mean(mean(decoder_tru==decoder_guess,3),1);
plot(t_mids,m,'Color',[0 .7 0]);
m_shuf = mean(mean(decoder_tru==decoder_shuf,3),1);
plot(t_mids,m_shuf,'Color',[.5 .5 .5]);

plot([0 0],[0 1],'k')
xlim([-750,1500])
ylim([.2 .6])
xlabel('time from pics (ms)')
ylabel('value decoder accuracy')

title({session,[region,' (PC exp var ',num2str(PCthresh),')']})

%% show conf mat at peak accuracy
[~,peak] = max(m);
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

