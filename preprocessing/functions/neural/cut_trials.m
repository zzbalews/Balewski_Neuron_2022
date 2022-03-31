function [sync_cut] = cut_trials(datavector,sync_entries)
    
datavector(end+1) = NaN;

% point incomplete trials to NaN;
sync_entries(isnan(sync_entries)) = length(datavector);

% cut trials
sync_cut = datavector(sync_entries);


end

