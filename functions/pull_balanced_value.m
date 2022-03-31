function [trialsLOO, setsize, N] = pull_balanced_value(minidata, keep_tr, niters)

% trial info
forced_tr = find(keep_tr);

lever = minidata.lever(keep_tr);
jpg = minidata.jpg(keep_tr,:);
lever_opts = [-1,1];

% all forced options: 1-16 on each side
forced_opts = [(1:16)', 17*ones(16,1); 17*ones(16,1), (1:16)'];
setsize = size(forced_opts,1); % # unique forced trials

% identify each trial as particular trial type
[~,idx] = ismember(jpg,forced_opts,'rows');

% find max number of balanced sets
counts = zeros(setsize,1);
for u = 1:setsize
    counts(u) = sum(idx==u);
end

N = min(counts);

% build balanced sets
trialsLOO = cell(niters,2); % col1 = balanced set; col2 = unused trials

for i = 1:niters
   
    % make nice sets
    trialsLOO{i,1} = zeros(setsize,N);
    for u = 1:setsize
    temp = find(idx==u);
    selected = forced_tr(randsample(temp,N));
    
    trialsLOO{i,1}(u,:) = selected;
    end
    
    % find unused trials
    used = reshape(trialsLOO{i,1},1,[]);
    unselected = forced_tr(~ismember(forced_tr,used));
    
    trialsLOO{i,2} = unselected';

end

