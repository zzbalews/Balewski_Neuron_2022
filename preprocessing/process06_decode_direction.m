%% train LDA decoder (direction)
% all preprocessing scripts run for this recording session: George rec18 02-11-2021
% swap out session details and file paths to run for other sessions


%% setup paths

addpath ../functions/
addpath functions/
addpath functions/decoding/
addpath functions/unit/

bhvdir = '../data_processed'; % where to find nice bhv data
pl2dir = '../data_processed'; % where to find smoothed/cut firing rates

% edit these for specific session
subject = 'George';
experiment = 'decoding_CN+OFC';
session_date = '02-11-2021';
session_num = 'George00_rec18';

dir_out = '../output_decoding'; % where to save decoding outputs

% which region?
region = 'CN'; % dir decoding doesn't work for OFC; just run this on CN data

% which trial period?
trialperiod = 'pics'; % can be 'pics' or 'choice'

which_sources = {'SPKfrnorm_units','LFP01delta'};


%% get filenames
session_name = strjoin({session_num,strrep(session_date,'-','')},'_');

% bhv
bhvfile = 'amntprob4x4_George_2021-02-11_clean.mat';

% spk -- smoothed
temp = dir([pl2dir,'/',session_name,'*w20_s5_',trialperiod,'*',region,'*']);
spk20 = {temp.name};


%% extract useful bhv variables
minidata = get_useful_bhv(bhvdir,bhvfile);

%% do LOO value decoding on free trials
decoder = 'direction';
niters = 25;

decode_output = do_free_LOO_fixPCA(minidata,pl2dir,spk20,dir_out,decoder,niters,trialperiod,which_sources);

%% show post prob/states
show_free_dir_decoding(decode_output,minidata,decoder)

