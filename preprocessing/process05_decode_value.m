%% train LDA decoder (value)

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
region = 'OFC'; % can be 'OFC' or 'CN'

% which trial period?
trialperiod = 'pics'; % can be 'pics' or 'choice'

%% get filenames
session_name = strjoin({session_num,strrep(session_date,'-','')},'_');

% bhv
bhvfile = 'amntprob4x4_George_2021-02-11_clean.mat';

% spk -- smoothed
temp = dir([pl2dir,'/',session_name,'*w200_s50_',trialperiod,'*',region,'*']);
spk200 = {temp.name};

temp = dir([pl2dir,'/',session_name,'*w20_s5_',trialperiod,'*',region,'*']);
spk20 = {temp.name};

temp = dir([pl2dir,'/',session_name,'*slice_w400*',trialperiod,'*',region,'*']);
slice400 = {temp.name};


%% extract useful bhv variables
minidata = get_useful_bhv(bhvdir,bhvfile);

%% do LOO value decoding on forced trials
decoder = 'value';
niters = 20;

% train/test on forced; LOO accuracy to find peak
do_forced_LOO_fixPCA(minidata,pl2dir,spk200,dir_out,decoder,niters,trialperiod);

%% train on forced slice; apply to free trials
decoder = 'value';
niters = 20;

% train on forced, apply weights to free
do_free_apply(minidata,pl2dir,spk20,slice400,dir_out,decoder,niters);

%% show training accuracy on 100-500 ms pics
decoder = 'value';
niters = 20;

% do LOO accuracy again, just on the same slice used for free trials
get_value_LOOpics(minidata,pl2dir,slice400,dir_out,decoder,niters);

