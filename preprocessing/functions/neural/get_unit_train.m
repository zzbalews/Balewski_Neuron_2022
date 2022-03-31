function  [spiketrain_smooth, unit_names, raw_fr] = get_unit_train(...
    do_zscore, ts, unit_names, chname, u, maxT)

% round spike times to nearest ms
ts_ms = round(ts*1000);
ts_ms = ts_ms(ts_ms>0); % remove if spike in 0th ms by rounding

% convert timestamps to spike train
spiketrain = ts_to_train(ts_ms, maxT);

% boxcar smooth
% spiketrain_smooth = 1000 * fastsmooth(spiketrain, 49, 1, 1);
spiketrain_smooth = spiketrain;

% normalized, if needed
if do_zscore
    mu = mean(spiketrain_smooth);
    sig = std(spiketrain_smooth);
    spiketrain_smooth = (spiketrain_smooth-mu)/sig;
end

% create unique unit label
unit_names = add_unit_label(unit_names, chname, u);

% raw firing rate
raw_fr = sum(spiketrain)/(length(spiketrain)/1000); %Hz

end

