%% single unit firing rates --> do glms
% all preprocessing scripts run for this recording session: George rec18 02-11-2021
% swap out session details and file paths to run for other sessions

%% setup paths

addpath ../functions/
addpath functions/
addpath functions/unit/

bhvdir = '../data_processed'; % where to find nice bhv data
pl2dir = '../data_processed'; % where to find smoothed/cut firing rates
pl2dir_temp = '../data_temp'; % where to find midprocessing neural data

dir_out = '../output_regressions'; % where to save regression outputs

% edit these for specific session
subject = 'George';
experiment = 'decoding_ACC+OFC';
session_date = '02-11-2021';
session_num = 'George00_rec18';

which_val_model = 'expval_bin';

% change 'trialperiod' and 'sigrange' if want to look at sig neurons before choice
trialperiod = 'pics'; %'pics' or 'choice'
sigrange = [0,500]; %ms from sync event: [0,500] or [-500, 0]

%% get filenames
session_name = strjoin({session_num,strrep(session_date,'-','')},'_');

% bhv
bhvfile = 'amntprob4x4_George_2021-02-11_clean.mat';

% spk -- smoothed
temp = dir([pl2dir,'/',session_name,'*w100_s25_',trialperiod,'*']);
spkfiles = {temp.name};

% spk -- raw (just to get firing rate for the whole session)
temp = dir([pl2dir_temp,'/',session_name,'*',trialperiod,'-SPKfr.mat']);
spkraw = temp(1).name;
units = get_unit_fr(pl2dir_temp,spkraw);

%% save glms

% extract useful variables from bhv
minidata = get_useful_bhv(bhvdir,bhvfile);

% do glm on all spk files, illustrate
for f = 1:length(spkfiles)
    
    % all trials: firing rate ~ trial type + lever dir + max val (+trial #)
    glmfile = do_glms(pl2dir, spkfiles{f}, minidata, units, dir_out, which_val_model, 'alltrials');
    show_unit_profiles({glmfile},'alltrials',sigrange); % # sig units & overall firing rates
    
end
