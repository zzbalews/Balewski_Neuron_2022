%% show behavior results
% load all data
% empirical behavior = choice, RT, saccades
% model results = psychometric fit


%% paths
addpath functions_bhv/
addpath functions_modelfit/
addpath ../functions/

fname_sessions = '../info/recording_sessions.csv';

dir_data = '/media/StorageHDD/WallisLabData';
dir_bhv = 'bhv';
dir_pl2 = 'pl2_cut';

dir_stims = '../info';

%% get session info
ses_info = readtable(fname_sessions,'Format','%s%s%s%s%d%d');

%% load all bhv data
bhvdata = load_cleanbhv(ses_info, dir_data, dir_bhv);

%% ~~~~ summary stats over raw behavior (empirical) ~~~~

%% summary stats for # sessions, # trials,s etc. (free, forced separate)
liststats_empirical_bhv(bhvdata);

%% heatmap: show images by subject & label with valbin
% Fig. 1b
show_16pics_jpg(bhvdata,dir_stims,'jpg');

%% heatmap: how often is each picture chosen? (free trials only)
% Fig. 1c
show_16pics_heatmap(bhvdata,'choice','free');

%% histogram: show overall RT, 1st sacc time, 2nd sacc time distributions (fre & forced)
% Fig. 1d
show_responsetime_hists(bhvdata);
set(gcf,'Position',[100 100 400 225])

%% ~~~~ summary stats over model fit (model) ~~~~

%% show psychometric fit & summary stats (% expvar, % chose high>low, etc)
% Fig. 1e
show_psychometric(bhvdata)
set(gcf,'Position',[100 100 500 225])
%% show lever RT vs max val (& vs delta val)
% Fig. 1g
show_responsetime_resid_vsvalue(bhvdata,'RT','valbin');
% show_responsetime_vsvalue(bhvdata,'RT','valbin');
set(gcf,'Position',[100 100 475 300])
%% show 1st sacc time vs max val (& vs delta val)
% Fig. 1h
show_responsetime_resid_vsvalue(bhvdata,'sacctime','valbin');
% show_responsetime_vsvalue(bhvdata,'sacctime','valbin');
set(gcf,'Position',[100 100 475 300])

%% value bin breakdown
% Fig. 1f
show_16pics_jpg(bhvdata,dir_stims,'valbin');


