function [lfp_range,lfp_mag] = loadlfp(fname,which_period)

if strcmp(which_period,'pics')
    load(fname,'pics_range','pics_mag')
    lfp_range = pics_range;
    lfp_mag = pics_mag;
elseif strcmp(which_period,'choice')
    load(fname,'choice_range','choice_mag')
     lfp_range = choice_range;
    lfp_mag = choice_mag;
elseif strcmp(which_period,'sacc1')
    load(fname,'sacc1_range','sacc1_mag')
     lfp_range = sacc1_range;
    lfp_mag = sacc1_mag;
elseif strcmp(which_period,'sacc2')
    load(fname,'sacc2_range','sacc2_mag')
     lfp_range = sacc2_range;
    lfp_mag = sacc2_mag;
end

end

