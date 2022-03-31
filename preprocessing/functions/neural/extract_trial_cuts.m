function [pics_entries, choice_entries, sacc1_entries, sacc2_entries] = extract_trial_cuts(...
    pics_range,choice_range,sacc1_range,sacc2_range,bhv,nan_placeholder)
% from bhv, get timestamps for trials sync to pics and choice
% load bhv file
load(bhv);

% number trials
ntr = length(trialinfo.TrialNumber);

% safety check
if isempty(pics_range)
    pics_range = [0,0];
end
if isempty(choice_range)
    choice_range = [0,0];
end
if isempty(sacc1_range)
    sacc1_range = [0,0];
end
if isempty(sacc2_range)
    sacc2_range = [0,0];
end

% time within trial
pics_idx = pics_range(1):pics_range(2);
choice_idx = choice_range(1):choice_range(2);
sacc1_idx = sacc1_range(1):sacc1_range(2);
sacc2_idx = sacc2_range(1):sacc2_range(2);

% number of time points
picsN = length(pics_idx);
choiceN = length(choice_idx);
sacc1N = length(sacc1_idx);
sacc2N = length(sacc2_idx);

% time across trials: @pics or @choice
pics_t = round(trialinfo.PL2picson);%round(pl2_info(:,5)); %ms
choice_t = round(trialinfo.PL2choice);%round(pl2_info(:,6)); %ms
sacc1_t = round(trialinfo.PL2sacc(:,1));
sacc2_t = round(trialinfo.PL2sacc(:,2));

% make matrix with all timepoints to pull: rows = trials, columns = time w/in trial
pics_entries = repmat(pics_t,1,picsN) + repmat(pics_idx,ntr,1);
choice_entries = repmat(choice_t,1,choiceN) + repmat(choice_idx,ntr,1);
sacc1_entries = repmat(sacc1_t,1,sacc1N) + repmat(sacc1_idx,ntr,1);
sacc2_entries = repmat(sacc2_t,1,sacc2N) + repmat(sacc2_idx,ntr,1);

% set incomplete trials as 1+length of filtered LFP signal (will be set to NaN)
pics_entries(pics_t==0,:) = nan_placeholder;
choice_entries(choice_t==0,:) = nan_placeholder;
sacc1_entries(sacc1_t==0,:) = nan_placeholder;
sacc2_entries(sacc2_t==0,:) = nan_placeholder;

end

