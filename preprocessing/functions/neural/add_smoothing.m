function smth_sets = add_smoothing(smth_sets,smth_type,w,s,mid,which_period,regions)

% check inputs are sensible
if strcmp(smth_type,'boxcar')
    if isempty(w) | isempty(s) | ~isempty(mid)
        error('add_smoothing error: boxcar inputs are not sensible')
    end
    
elseif strcmp(smth_type,'slice')
    if isempty(w) | ~isempty(s) | isempty(mid)
        error('add_smoothing error: slice inputs are not sensible')
    end
    
else
    error('add_smoothing error: smooth type not recognized')
end


% add to smoothing list
if length(fieldnames(smth_sets))==0
    n = 1;
else
    n = length(smth_sets)+1;
end

smth_sets(n).type = smth_type;
smth_sets(n).w = w;
smth_sets(n).s = s;
smth_sets(n).mid = mid;

if strcmp(which_period,'both')
    trialperiods = {'pics','choice'};
elseif strcmp(which_period,'pics') | strcmp(which_period,'choice') | ...
        strcmp(which_period,'sacc1') | strcmp(which_period,'sacc2')
    
    trialperiods = {which_period};
else
    error('did not recognize trial period: pics, choice, both')
end
smth_sets(n).trialperiods = reshape(trialperiods,1,[]);

smth_sets(n).regions = reshape(regions,1,[]);

% add expected filename fragments
if strcmp(smth_type,'boxcar')
    smth_sets(n).label = ['_smoothed_',smth_type,'_w',num2str(w),'_s',num2str(s)];
elseif strcmp(smth_type,'slice')
    smth_sets(n).label = ['_smoothed_',smth_type,'_w',num2str(w),'_T',num2str(mid)];
end


end

