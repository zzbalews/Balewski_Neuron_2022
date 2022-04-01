% compare CN direction decoding by first OFC value state
%% paths
addpath functions/
addpath ../figures_behavior/functions_bhv/
addpath ../functions/

fname_sessions = '../info/recording_sessions.csv';

dir_data = '/media/StorageHDD/WallisLabData';
dir_bhv = 'bhv';

% select whether looking at OFC or CN value decoding

dir_decode_val = 'decoding_SAVEvalueOFC';
dir_decode_dir = 'decoding_SAVEdirectionCN';

thresh = 0.001;

dirssd = '/media/StorageSSD/Dropbox/Documents/WallisLab/analysis/data_working';


%% get session info

ses_info = readtable(fname_sessions,'Format','%s%s%s%s%d%d');

ses_info = ses_info(ses_info.OFCvaluedecoding==1,:);

%% load bhv data
bhvdata = load_cleanbhv(ses_info, dir_data, dir_bhv);

%% load value decoder
valdata = load_valdecoder(ses_info,bhvdata,dirssd, dir_decode_val, 0);

%% load direction deocoder
dirdata = load_dirdecoder('pics', ses_info, dirssd, dir_decode_dir, bhvdata, 0);

%% get postprob dist from fixation period
[dirdata] = get_fix_dist(bhvdata, dirdata, []);

%% restrict to trials: ch > unch value bin
keep_tr = bhvdata.chose_high>0 & ...
    bhvdata.lever~=0 & ...
    bhvdata.trialtype==2;% & ismember(bhvdata.session,deltasessions);

thresh = 0.0001;

%% split trials by first OFC value state
% Fig. 6e

[~, t0] = min(abs(valdata.t_mids-0));
[~, t_end] = min(abs(valdata.t_mids-800));

ntr = length(bhvdata.TrialNumber);

first_OFC_state = zeros(ntr,1);

for tr = find(keep_tr)'
    
    idx_ch = find(valdata.ch_state(tr, t0:t_end), 1);
    idx_unch = find(valdata.unch_state(tr, t0:t_end), 1);
   
    if isempty(idx_ch) & isempty(idx_unch)
        first_OFC_state(tr) = 0;
        
    elseif ~isempty(idx_ch) & isempty(idx_unch)
        first_OFC_state(tr) = 1;
        
    elseif isempty(idx_ch) & ~isempty(idx_unch)
        first_OFC_state(tr) = -1;
        
    elseif ~isempty(idx_ch) & ~isempty(idx_unch)
        
        [~, temp] = min([idx_ch, idx_unch]);
        
        if temp==1
            first_OFC_state(tr) = 1;
        elseif temp==2
            first_OFC_state(tr) = -1;
        end
    end
    
end

bhvdata.first_OFC = first_OFC_state;
tr_sampled = sample_matched_trials(bhvdata,keep_tr,1,0,'which1stOFCval');
show_dirdec_byfirstgaze_mini(bhvdata,dirdata,tr_sampled,'direction',keep_tr);
set(gcf,'Position',[50 50 560 250])
