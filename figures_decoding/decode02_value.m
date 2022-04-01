% show OFC and CN value decoding results

%% setup
addpath functions/
addpath ../figures_behavior/functions_bhv/
addpath ../functions/

fname_sessions = '../info/recording_sessions.csv';

dir_data = '/media/StorageHDD/WallisLabData';
dir_bhv = 'bhv';

% select whether looking at OFC or CN value decoding
region = 'OFC';
if strcmp(region,'OFC')
    dir_decode = 'decoding_SAVEvalueOFC';
    dir_decode_choice = 'decoding_SAVEvalueOFC_choice';
else
    dir_decode = 'decoding_SAVEvalueCN';
    dir_decode_choice = 'decoding_SAVEvalueCN_choice';
    
end

thresh = 0.001;

dirssd = '/media/StorageSSD/Dropbox/Documents/WallisLab/analysis/data_working';


%% get session info

ses_info = readtable(fname_sessions,'Format','%s%s%s%s%d%d');

if strcmp(region,'OFC')
    ses_info = ses_info(ses_info.OFCvaluedecoding==1,:);
else
    ses_info = ses_info(ses_info.CNvaluedecoding==1,:);
end

%% load bhv data
bhvdata = load_cleanbhv(ses_info, dir_data, dir_bhv);

%% load value decoder
get_LOO = 1;
[valdata,LOO] = load_valdecoder(ses_info, bhvdata, dirssd, dir_decode, get_LOO);
if strcmp(region,'OFC')
    [valdata_choice, ~] = load_valdecoder(ses_info, bhvdata, dirssd, dir_decode_choice, 0);
end

if get_LOO
    show_LOO_acc(ses_info,LOO)
end
%% restrict to trials ch > unch valuebin

keep_tr = bhvdata.chose_high>0 & ...
    bhvdata.lever~=0 & ...
    bhvdata.trialtype==2;

thresh = 0.001;
%% show avg ppd ch v unch v na
% Fig. S2
show_ppd(bhvdata,valdata,keep_tr,'pics');
set(gcf,'Position',[50 50 560 250])

show_ppd(bhvdata,valdata_choice,keep_tr,'choice');
set(gcf,'Position',[50 50 560 250])

%% log10(RT) ~ 1 + logit(ch ppd) + logit(unch ppd)
% Fig. 6c, 7c
show_glmRT_ppd(bhvdata,valdata,keep_tr);
set(gcf,'Position',[50 50 560 250])

%% count # ch/unch/na states & flips
% Fig. 6a-b, 7a-b
decoder = 'value';
count_states(bhvdata,valdata,keep_tr,decoder);
set(gcf,'Position',[50 50 560 250])

%% check if post prob different when 1st loc = lever vs not (matched trials)
% Fig. 6d, 7d
tr_sampled = sample_matched_trials(bhvdata,keep_tr,100,0,'whichsacc');
% set(gcf,'Position',[50 50 560 125])

show_dirdec_byfirstgaze(bhvdata,valdata,tr_sampled,'value','chosen',keep_tr);
set(gcf,'Position',[50 50 560 250])

show_dirdec_byfirstgaze(bhvdata,valdata,tr_sampled,'value','unchosen',keep_tr);
set(gcf,'Position',[50 50 560 250])
