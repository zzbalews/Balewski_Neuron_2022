function spiketrain = ts_to_train(ts_ms, maxT)
% convert spike timestamps to spike train

last_point = 1000 + max([maxT; ts_ms]);

spiketrain = zeros(1,last_point);

spiketrain(ts_ms) = 1;

end

