%% turn spike times into firing rates; bandpass LFP
% all preprocessing scripts run for this recording session: George rec18 02-11-2021
% swap out session details and file paths to run for other sessions

%% setup -- edit variables for session, then run!

addpath ../functions/
addpath functions 
addpath functions/neural/ % specific bhv
addpath functions/'Matlab Offline Files SDK'/ % plexon

% get location for each channel in recordings
channelmap = '../info/George_channel_source.csv'; % recording location for each channel

bhvdir = '../data_processed'; % where to find processed behavior .mat file
pl2dir = '../data_raw'; % where to find LFP and spike .pl2 files
pl2_temp = '../data_temp'; % temp dir to save neural data before aligning to trial events
pl2_out = '../data_processed'; % where to spit out processed spike and bandpassed LFP magnitude, smoothed and aligned to trial events

subject = 'George';
experiment = 'decoding_CN+OFC';
session_date = '02-11-2021';
session_num = 'George00_rec18';

% align to trial events: cut times
pics_range = [-2000,2000]; % ms range around picture onset
choice_range = [-2000,2000]; % ms range around choice

sync_opts = {'pics','choice'}; % name

%% get filenames
session_name = strjoin({session_num,strrep(session_date,'-','')},'_');

% bhv
bhvfile = 'amntprob4x4_George_2021-02-11_clean.mat';

% neural raw
spkfiles = {'George00_rec18_02112021-units.pl2'}; % list all spike files (*unit) for this session
lfpfile = 'George00_rec18_02112021-LFP.pl2';

%% get region information

channelID = get_channel_regions(channelmap,session_date,session_num);
regions = fieldnames(channelID); % by default, assuming do smoothing in both regions

restricted_channels = [];
for rr = 1:length(regions)
    restricted_channels = cat(1,restricted_channels,...
        table2array(channelID.(regions{rr})));
end

%% turn spike times into firing rates
for do_zscore = 0:1 % 0=raw, 1=zscore magnitudes
    for s = 1:2 % 1=sync to pictures on, 2=sync to choice
        which_sync = sync_opts{s};
        
        tic
        channels_with_spk = make_spk_cut(bhvdir,bhvfile,pl2dir,spkfiles,pl2_temp,...
            pics_range,choice_range,[],[],session_name,do_zscore,which_sync,...
            restricted_channels);
        toc
        
    end
end

%% filter LFP (common bands)
make_lfp_cut(bhvdir,bhvfile,pl2dir,lfpfile,pl2_temp,...
    pics_range,choice_range,[],[],session_name,channels_with_spk);

%% combine spk & lfp, split by region, do smoothing
% specify which boxcar smoothing @ each trial sync point

smth_sets = struct(); % collect smoothings to do

% boxcar: window 200ms stepped by 50ms; @ pics and choice
smth_sets = add_smoothing(smth_sets,'boxcar',200,50,[],'both',regions); 
% boxcar: window 100ms stepped by 25ms; @ pics and choice
smth_sets = add_smoothing(smth_sets,'boxcar',100,25,[],'both',regions);
% boxcar: window 20ms stepped by 5ms; @ pics and choice
smth_sets = add_smoothing(smth_sets,'boxcar',20,5,[],'both',regions);
% average in slice: 400ms window around 300ms after pics
smth_sets = add_smoothing(smth_sets,'slice',400,[],300,'pics',regions);
% average in slice: 400ms window around -200ms before choice
smth_sets = add_smoothing(smth_sets,'slice',400,[],-200,'choice',regions);

% will check if smoothed files already exist before repeating:
smooth_brain_data(pl2_temp, pl2_out, channelID, smth_sets, session_name)

