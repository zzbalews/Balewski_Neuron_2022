  %% show single unit regression results



%% paths
addpath functions/
addpath functions/venn/
addpath ../functions/
addpath ../figures_behavior/functions_bhv/

fname_sessions = '../info/recording_sessions.csv';

dir_data = '/media/StorageHDD/WallisLabData';
dir_bhv = 'bhv';
dir_pl2 = 'pl2_cut';
dir_out = '.';

experiment = 'decoding_CN+OFC';
which_value = 'pics_expval_bin'; 
dir_unitglms = ['singleunits/',which_value];
which_model = 'alltrials'; % 'alltrials' or 'free'

dirssd = '/media/StorageSSD/Dropbox/Documents/WallisLab/analysis/data_working';

%% get session info
ses_info = readtable(fname_sessions,'Format','%s%s%s%s%d%d');

%% load all unit glms
glmdata = load_glmresults(ses_info,dirssd, ...
    dir_unitglms, which_model);

%% display unit & LFP channel counts

%% show venn diagram of sig units (@ different sig levels?)
% Fig. 2b, 5b
out_prefix = [which_value,'_',which_model];
unitcounts = show_unit_venn(glmdata,[0 500],out_prefix);

%% show proportion sig units between regions
% Fig. 5c
show_unit_proportions(unitcounts);
set(gcf,'Position',[50 50 560 200])

%% show CPD beta heatmaps
% Fig. 2d-e, 5e
show_cpd_heatmaps(glmdata,out_prefix)

%% show example PSTHs
% Fig. 2c, 5d
%load bhv data
bhvdata = load_cleanbhv(ses_info, dir_data, dir_bhv);

example_units = [633,1,134,104,1546]; % CN & OFC examples
make_nice_psths(bhvdata,glmdata,example_units,dirssd,dir_pl2);

%% compare 1st sig bin (value units only) between CN & OFC
compare_firstsig(glmdata)

