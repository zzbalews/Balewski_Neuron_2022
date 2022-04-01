function get_value_LOOpics(minidata,pl2dir,spkslice,outdir,decoder,niters)

%% data details

temp = strsplit(spkslice{1},'_');
session = strjoin(temp(1:3),' ');
if length(spkslice)==1
    region = strrep(temp{end},'.mat','');
else
    region = strrep(temp{end}(2:end),'.mat','');
end


% outprefix_old = strjoin({[outdir,'/old'],[strrep(session,' ','_'),'_',region,'_value_free_w400']},'/');
outprefix = strjoin({outdir,[strrep(session,' ','_'),'_',region,'_value_LOOacc_w400']},'/');

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


%% load brain data
varnames = {'SPKfrnorm_units',...
    'LFP01delta','LFP02theta','LFP03alpha',...
    'LFP04beta','LFP05gamma','LFP06highgamma'};

% varnames = {'SPKfrnorm_units'};

brain_slice = cell(1,length(varnames));

unit_names = {};
for i = 1:length(spkslice)
    
    % load slice (forced; train)
    load(strjoin({pl2dir,spkslice{i}},'/'),...
        'SPKfrnorm_names','t_mids');
    unit_names = cat(1,unit_names,SPKfrnorm_names);
    
    for v = 1:length(varnames)
        slice = load(strjoin({pl2dir,spkslice{i}},'/'), varnames{v});
        brain_slice{v} = cat(3, brain_slice{v}, slice.(varnames{v}));

    end
end

nmids = length(t_mids);

%% cycle through iters, average 
PCthresh = 95;

track_label = nan(N,setsize,niters);
track_guess = nan(N,setsize,niters);
tic
for i = 1:niters
   
    % get forced trials in this iter
    
    for S = 1:N % cycle through LOO sets
        
        
        LOOset = trialsLOO{i,1};
        
       test_tr = reshape(LOOset(:,S),[],1);
       input_tr = reshape(LOOset(:,[1:(S-1),(S+1):end]),[],1);
       
       test_label = max(truelabel(test_tr,:),[],2);
       input_label = max(truelabel(input_tr,:),[],2);
    
    % get top pcs for slice
    input_pcs = [];
    test_pcs = [];
    for freq = 1:length(varnames)
        snippet_train = permute(brain_slice{freq}(input_tr,:,:),[1,3,2]);
        [coeff, score, ~, ~, expvar, mu] = pca(snippet_train);
        PC95 = find(cumsum(expvar)>PCthresh,1);
        
        input_pcs = cat(2, input_pcs, score(:,1:PC95));

    % use same pc projections for smooth
    snippet_apply = permute(brain_slice{freq}(test_tr,:,:),[1,3,2]);
    
    term1 = snippet_apply - repmat(mu,setsize,1,nmids);
    term2 = repmat(coeff(:,1:PC95),1,1,nmids);
    
% % % %     apply_pcs = nan(ntr,PC95,nmids);
% % % %     for m = 1:nmids
% % % %         apply_pcs(:,:,m) = term1(:,:,m) * term2(:,:,m);
% % % %     end
    test_pcs = cat(2, test_pcs, pagemtimes(term1,term2)); % ntr x pcs x nmids
    end
    
    % train decoder
    mdl_iter = fitcdiscr(input_pcs,input_label); % use default LDA
    [test_guess] = predict(mdl_iter,test_pcs);
    
    
    track_label(S,:,i) = test_label;
    track_guess(S,:,i) = test_guess;
    
    end
    
end
toc
%% average LOO accuracy
m = mean(track_label==track_guess,3);
mean(reshape(m,[],1));

%% conf mat
C = confusionmat(reshape(track_label,[],1),...
    reshape(track_guess,[],1),'order',levels);
C = C';

total = repmat(sum(C,1),nlevels,1);
conf_mat = C./total;

figure;
imagesc(conf_mat)
daspect([1 1 1])

colormap(hot)
colorbar
caxis([.25 .8])

title([strrep(session,'_','\_'),': ',num2str(round(mean(reshape(m,[],1)),2))])
%% save new files

print([outprefix,'.png'],'-dpng')


 save([outprefix,'.mat'],'conf_mat','levels','track_label','track_guess');

end
