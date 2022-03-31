function [trialsLOO_free, setsize, N, leftover] = pull_balanced_direction(minidata, truelabel, keep_tr, niters)
% make balanced sets of trials using value and direction for free only
% (use ch>=unch only for training; track all unused trials to apply
% decoder)

% trial info
free_tr = find(keep_tr);

val = minidata.valbin;
val(isnan(val)) = 0;

lever_opts = [-1,1];

%% decide unique trial types for balancing

% free: only ch > unch on L & R (can also try ch>=unch by changing line 20)
free_trialtypes = [];
for v = 1:4
     for u = (v+1):4 %for u = v:4 %
        temp = [v, u, 1];
        free_trialtypes(end+1,:) = temp;
    end
end

temp = [free_trialtypes(:,[2,1]), -free_trialtypes(:,3)];
free_trialtypes = cat(1, free_trialtypes, temp);

nfreetypes = size(free_trialtypes,1);

setsize = nfreetypes;

%% identify each trial as one of these trial types (or 0 = do not use)

[~, free_idx] = ismember([val(free_tr,:),truelabel(free_tr)],...
    free_trialtypes, 'rows');

%% find larged N (# of balanced sets)

counts_free = zeros(nfreetypes,1);
for u = 1:nfreetypes
    counts_free(u) = sum(free_idx==u);
end

N = min(counts_free);

%% build sets of trials

trialsLOO_free = cell(niters,2); % col1 = balanced sets; col2 = unused trials
for i = 1:niters
    
    % select balanced trials for LOO
    trialsLOO_free{i,1} = zeros(setsize, N);
    for u = 1:nfreetypes
        temp = find(free_idx==u);
        selected = free_tr(randsample(temp, N));
        
        trialsLOO_free{i,1}(u,:) = selected;
    end
    
    % save remaining free trials, not used for LOO
    remaining = free_tr(~ismember(free_tr, reshape(trialsLOO_free{i},[],1)));
    ntr = length(remaining);
    bonus = ceil(ntr/N)*N-ntr;
    
    trialsLOO_free{i,2} = reshape([remaining(randperm(ntr)); repmat(NaN,bonus,1)],[],N);
end

leftover = size(trialsLOO_free{i,2},1);


end