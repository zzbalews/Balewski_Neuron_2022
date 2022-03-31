function [timestamp, t_mids] = get_midpoints(smoothtype, ...
    slice, t_range, w, s)

% timestamp of columns in neural data
timestamp = t_range(1):t_range(2);

if strcmp(smoothtype,'slice')
    t_mids = slice(1);
    
elseif strcmp(smoothtype,'boxcar')
    mids = floor(w/2):s:(length(timestamp)-floor(w/2));
    t_mids = timestamp(mids);
    
end

end

