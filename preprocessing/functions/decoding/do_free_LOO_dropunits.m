function do_free_LOO_dropunits(minidata, pl2dir, spksmooth, outdir, ...
    decoder, niters, inc, trange)

%% data details

temp = strsplit(spksmooth{1},'_');
session = strjoin(temp(1:3),' ');
region = strrep(temp{end}(2:end),'.mat','');

outfile_prefix = strjoin({outdir,[strrep(session,' ','_'),'_',region,'_direction_neurondrop']},'/');

%% set up decoder details
if strcmp(decoder,'direction')
    
    levels = [-1,1];
    nlevels = 2;
    
    truelabel = minidata.lever;

else
    disp('need to set up this decoder! STOP')
    return
end


%% select trials

keep_tr = minidata.trialtype==2 & minidata.lever~=0;

% make balanced sets
[trials_balanced, setsize, N, leftover] = pull_balanced_direction(minidata, keep_tr, niters);

trialsLOO = trials_balanced(:,1); % balanced sets to use for LOO

ntr = setsize * N;


%% load brain data
braindata = [];
for f = 1:length(spksmooth)
    load(strjoin({pl2dir,spksmooth{f}},'/'),...
    't_mids','SPKfrnorm_units');
    
    braindata = cat(3,braindata,SPKfrnorm_units);
end

t_keep = find(t_mids>=trange(1) & t_mids<=trange(2));
t_keep = t_keep(1:4:length(t_keep));

% restrict to peak decoding period
braindata = braindata(:,t_keep,:);
t_mids = t_mids(t_keep);

nmids = length(t_mids);

nunits = size(braindata,3);

%% determine neuron sets
ncounts = inc:inc:nunits;

sampled_units = cell(1,length(ncounts));
nsampl = 50;
for i = 1:length(ncounts)
   
    temp = nan(nsampl,ncounts(i));
    for iter = 1:nsampl
        temp(iter,:) = randsample(nunits,ncounts(i));
    end
    
    sampled_units{i} = temp;
    
end

%% do decoding
best_acc = nan(length(ncounts),nsampl,nmids);

parfor C = 1:length(ncounts) % inc # units
for P = 1:nsampl % sample within inc
    tic
    current_units = sampled_units{C}(P,:);
    
    for m = 1:nmids % times
        [~,score,~,~,expvar] = pca(permute(braindata(:,m,current_units),[1,3,2]));
        PC95 = find(cumsum(expvar)>95,1);
        braindata_pc = score(:,1:PC95);
        
        slice_truelabel = nan(ntr, 1, niters);
        slice_guess = nan(ntr, 1, niters);
        
        for i = 1:niters % shuffled trial iterations
            for s = 1:N % LOO sets
                
                tr_train = reshape(trialsLOO{i}(:,[1:s-1, s+1:end]),[],1);
                tr_test = trialsLOO{i}(:,s);
                
                % brain data for trials
                input_train = braindata_pc(tr_train,:);
                input_test = braindata_pc(tr_test,:);
                
                % correct label for trials
                label_train = truelabel(tr_train);
                label_true = truelabel(tr_test);
                
                % do LDA
                try
                [guess_test, ~, ~] = classify(input_test, input_train, ...
                    label_train,'linear','empirical'); % true LOO
                catch
                    guess_test = nan(size(label_true));
                end
                % update slice
                idx_test = (s-1)*setsize + (1:setsize)';
                slice_truelabel(idx_test,1,i) = label_true;
                slice_guess(idx_test,1,i) = guess_test;
            end
        end
        
        % update accuracy
        best_acc(C,P,m) = mean(nanmean(slice_truelabel==slice_guess,3));
            
    end
    
   toc
        
end

end

save([outfile_prefix,'.mat'],'best_acc','ncounts','t_mids','nsampl')

%% figure to illustrate
figure;
hold on
x = repmat(ncounts',1,nsampl);

scatter(reshape(x,1,[])+2,reshape(best_acc,1,[]),'filled','MarkerFaceAlpha',.1,'MarkerFaceColor',[1 0 0])
plot(ncounts,mean(best_acc,2),'k','LineWidth',2)

plot(ncounts([1,end]),[.5 .5]','k:')
plot(ncounts([1,end]),[.85 .85]','k:')
xlabel('# CN units')
ylabel('max direction accuracy')
ylim([.45 0.9])
xlim([0,ncounts(end)+diff(ncounts([1,2]))])

print([outfile_prefix,'.png'],'-dpng')
end

