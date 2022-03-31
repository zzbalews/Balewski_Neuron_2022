function minidata = get_useful_bhv(bhvdir,bhvfile)

fname = strjoin({bhvdir,bhvfile},'/');
load(fname)

minidata = struct();

vars = {'TrialNumber','lever','rt','ifreward','trialtype','jpg','saccloc'};
for v = 1:length(vars)
    minidata.(vars{v}) = trialinfo.(vars{v});
end

minidata.lever(isnan(minidata.lever)) = 0;

% use expected value model for subjective value!
minidata.subjval = trialinfo.subjval_expval;
minidata.valbin = trialinfo.valbin_expval;

% label chose high (by subjval)
[~,temp] = max(minidata.subjval,[],2);
best_loc = temp*2-3;
minidata.chosehigh = best_loc==minidata.lever;

% do 1st & 2nd saccade match lever?
minidata.sacc_to_lever = double(trialinfo.saccloc==trialinfo.lever);
minidata.sacc_to_lever(isnan(trialinfo.saccloc)) = NaN;


end

