function [w_fit_best,ll_fit_best] = best_fminsearch(loglikefunc,w_names)

options = struct();
options.MaxFunEvals = 1e5;
options.Display = 'off';

nparams = length(w_names);
w0 = [linspace(-1,1,3)',zeros(3,nparams-1)];

for i = 1:size(w0,1)
    [w_fit,ll_fit] = fminsearch(loglikefunc,w0(i,:),options);
    track.w_fit(i,:) = w_fit;
    track.ll_fit(i,:) = ll_fit;
end

[~,best] = min(track.ll_fit);

w_fit_best = track.w_fit(best,:);
ll_fit_best = track.ll_fit(best);

end

