function trial_timestamp = pull_eventcode(codes,event)


trial_timestamp = cellfun(@(x) x.CodeTimes(find(x.CodeNumbers==event,1,'first')),...
    codes,'UniformOutput',false);

idx_na = cellfun(@(x) isempty(x), trial_timestamp);

trial_timestamp(idx_na) = {NaN};

trial_timestamp = [trial_timestamp{:}]';




end

