function [timestamp, midpoints] = get_smth_points(smth_sets,missing_smooths,spkfr_range)

% get midpoints for all smoothing that needs to be done
midpoints = cell(size(missing_smooths));

        for i = 1:length(missing_smooths)
            iter = missing_smooths(i);
            
            [timestamp,t_mids] = get_midpoints(smth_sets(iter).type,...
                smth_sets(iter).mid,spkfr_range,...
                smth_sets(iter).w,smth_sets(iter).s);
            
            midpoints{i} = t_mids;
        end


end

