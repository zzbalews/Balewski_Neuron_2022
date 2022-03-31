function [chi2,df,chi2p] = do_chi2(observed)
% [chi2,df,chi2p] = do_chi2(observed)
% do chi2 test on 'observed' array


% dims of observed matrix
[nrows,ncols] = size(observed);
df = (nrows-1)*(ncols-1);

% calculate expected
expected = repmat(sum(observed),nrows,1) .* repmat(sum(observed,2),1,ncols) / ...
    sum(reshape(observed,1,[]));

% get chi2 stat
chi2 = sum(reshape(((observed-expected).^2)./expected,1,[]));

% get p value
chi2p = 1 - chi2cdf(chi2,df);

end

