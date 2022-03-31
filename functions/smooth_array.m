function X_smooth = smooth_array(X, t, t_mids, w)
% smooth X (trials x time x units) in time dimension at midpoints t_mids
% with window=w


% setup smoothed array
[ntr, ~, nunits] = size(X);
nmids = length(t_mids);
X_smooth = nan(ntr,nmids,nunits);

for m = 1:nmids
    
    if isnan(t_mids(m))
        continue
    end
    
    % get timepoints in window
    rng = t_mids(m) + [1-floor(w/2), floor(w/2)];
    idx = t>=rng(1) & t<=rng(2);
    
    % save mean
    X_smooth(:,m,:) = mean(X(:,idx,:),2);
end

end

