% compare CN direction decoding by first OFC value state
%% paths
addpath functions/
addpath ../figures_behavior/functions_bhv/
addpath ../functions/

fname_sessions = '../info/recording_sessions.csv';

dir_data = '/media/StorageHDD/WallisLabData';
dir_bhv = 'bhv';

% select whether looking at OFC or CN value decoding

dir_decode_valOFC = 'decoding_SAVEvalueOFC';
dir_decode_valCN = 'decoding_SAVEvalueCN';

thresh = 0.001;

dirssd = '/media/StorageSSD/Dropbox/Documents/WallisLab/analysis/data_working';


%% get session info

ses_info = readtable(fname_sessions,'Format','%s%s%s%s%d%d');

ses_info = ses_info(ses_info.OFCvaluedecoding==1 & ses_info.CNvaluedecoding==1,:);

%% load bhv data
bhvdata = load_cleanbhv(ses_info, dir_data, dir_bhv);

%% load value decoders
valdataOFC = load_valdecoder(ses_info,bhvdata,dirssd, dir_decode_valOFC, 0);
valdataCN = load_valdecoder(ses_info,bhvdata,dirssd, dir_decode_valCN, 0);

%% restrict to trials: ch > unch value bin
keep_tr = bhvdata.chose_high>0 & ...
    bhvdata.lever~=0 & ...
    bhvdata.trialtype==2;% & ismember(bhvdata.session,deltasessions);

thresh = 0.0001;

%% correlations between OFC and CN value decoding
% Fig. 7e

checkout_decoder_correlations(bhvdata,valdataOFC,valdataCN,keep_tr,'chosen');
set(gcf,'Position',[50 50 560 250])
checkout_decoder_correlations(bhvdata,valdataOFC,valdataCN,keep_tr,'unchosen');
set(gcf,'Position',[50 50 560 250])