function [spkfr_allunits, spkfr_range, unit_names] = loadspk(pl2processed_dir,which_spk,which_period)

%% find matching file for 'which_spk'

temp = dir(strjoin({pl2processed_dir,['*',which_period,'*',which_spk,'.mat']},'/'));

% should see only 1; give error if multiple
if isempty(temp)
    error([which_spk,' file missing']);
elseif length(temp)>1
    error([which_spk,' file conflict']);
end

SPKfr = strjoin({pl2processed_dir,temp.name},'/');

%% load variables corresponding to 'which_period'

if strcmp(which_period,'pics')
    load(SPKfr,'pics_allunits','pics_range','unit_names');
    spkfr_allunits = pics_allunits;
    spkfr_range = pics_range;
elseif strcmp(which_period,'choice')
    load(SPKfr,'choice_allunits','choice_range','unit_names');
    spkfr_allunits = choice_allunits;
    spkfr_range = choice_range;
elseif strcmp(which_period,'sacc1')
    load(SPKfr,'sacc1_allunits','sacc1_range','unit_names');
    spkfr_allunits = sacc1_allunits;
    spkfr_range = sacc1_range;
elseif strcmp(which_period,'sacc2')
    load(SPKfr,'sacc2_allunits','sacc2_range','unit_names');
    spkfr_allunits = sacc2_allunits;
    spkfr_range = sacc2_range;
end

end

