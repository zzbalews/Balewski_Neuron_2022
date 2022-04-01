% show CN direction decoding results

%% setup
addpath functions/
addpath ../figures_behavior/functions_bhv/
addpath ../functions/

fname_sessions = '../info/recording_sessions.csv';

dir_data = '/media/StorageHDD/WallisLabData';
dir_bhv = 'bhv';

dir_decode_pics = 'decoding_SAVEdirectionCN';
dir_decode_choice = 'decoding_SAVEdirectionCN_CHOICE';

thresh = 0.001;

dirssd = '/media/StorageSSD/Dropbox/Documents/WallisLab/analysis/data_working';
%% get session info
ses_info = readtable(fname_sessions,'Format','%s%s%s%s%d%d');

%% load bhv data
bhvdata = load_cleanbhv(ses_info, dir_data, dir_bhv);

%% load direction decoder
get_LOO = 1;
[dirdata_pics,LOO_pics] = load_dirdecoder('pics',ses_info, dirssd, dir_decode_pics, bhvdata, get_LOO);
[dirdata_choice,LOO_choice] = load_dirdecoder('choice',ses_info, dirssd, dir_decode_choice, bhvdata, get_LOO);

%% show LOO decoding accuracy
if get_LOO
    % Fig. 3a
show_dir_LOO_acc('pics',LOO_pics, bhvdata);

% Fig. 3b
compare_LOO_peaks(LOO_pics,LOO_choice)
end

%% get postprob dist from fixation period
[dirdata_pics,dirdata_choice] = get_fix_dist(bhvdata, dirdata_pics, dirdata_choice);

%% restrict trials to free & ch >= unch valuebin
keep_tr = bhvdata.chose_high>=0 & ...
    bhvdata.lever~=0 & ...
    bhvdata.trialtype==2;

%% show heatmat with chosen dir
% Fig. 3c
show_dirdec_heatmap('pics',bhvdata,dirdata_pics,keep_tr);
set(gcf,'Position',[50 50 560 1300])

%% show dir dec vs maxval (& GLM)
% Fig. 3d-e
show_dirdec_bymaxval(bhvdata,dirdata_pics,keep_tr,thresh,0);
set(gcf,'Position',[50 50 560 250])

%% show GLM RT ~ postprob (+maxval)
% Fig. 3f
show_dirdec_RTglm(bhvdata,dirdata_pics,keep_tr,thresh);
set(gcf,'Position',[50 50 560 250])

%% count # ch/unch/na states & flips
% Fig. 4d-e
decoder = 'direction';
count_states(bhvdata,dirdata_pics,keep_tr,decoder);
set(gcf,'Position',[50 50 560 250])

%% distribution of "state" onsets, using fix dist for cutoff
% Fig. 4a
show_dirdec_firststate(bhvdata,dirdata_pics,keep_tr);
set(gcf,'Position',[50 50 560 250])

%% check if post prob different when 1st loc = lever vs not (matched trials)
% Fig. 4c
tr_sameval = bhvdata.chose_high==0 & ...
    bhvdata.lever~=0 & ...
    bhvdata.trialtype==2;

tr_sampled = sample_matched_trials(bhvdata,tr_sameval,10000,0, 'whichsacc');
% set(gcf,'Position',[50 50 560 125])

show_dirdec_byfirstgaze(bhvdata,dirdata_pics,tr_sampled,'direction',[],tr_sameval);
set(gcf,'Position',[50 50 560 250])
