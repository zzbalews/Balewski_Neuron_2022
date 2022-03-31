%% convert .bhv or .bhv2 --> useful struct in .mat
% all preprocessing scripts run for this recording session: George rec18 02-11-2021
% swap out session details and file paths to run for other sessions

%% setup info

addpath ../functions/
addpath functions
addpath functions/behavior/ % specific bhv 
addpath functions/'Matlab Offline Files SDK'/ % plexon 

bhvdir = '../data_raw'; % where to find .bhv or .bhv2 files
pl2dir = '../data_raw'; % where to find LFP and spike .pl2 files
outdir = '../data_processed/'; % where to spit out processed files

% edit these for specific session (al
subject = 'George';
experiment = 'decoding_CN+OFC';
session_date = '02-11-2021';
session_num = 'George00_rec18';

% which psychometric model (linear/hyperbolic/exponential; + saccade)
psych_models = {'linear discount + first saccade', 'expected value + first saccade'};

%% get filenames
session_name = strjoin({session_num,strrep(session_date,'-','')},'_');

% bhv input
bhvfiles = {'amntprob4x4_George_2021-02-11.bhv2'}; % list all file names in order if multiple (had to restart task on some days)
if ~isempty(strfind(bhvfiles{1},'.bhv2'))
    % .bhv2 file, from ML2
    filetype = '.bhv2';
else
    % .bhv file, from ML1
    filetype = '.bhv';
end

% bhv output
temp = strsplit(bhvfiles{1},'_');
nicefile = strrep(strjoin([temp(1:3),{'clean.mat'}],'_'),filetype,'');

% lfp
lfpfile = 'George00_rec18_02112021-LFP.pl2';


%% pull useful trial features, add psychometric fits & save eye tracker data @ 1Hz

if strcmp(filetype,'.bhv') % ML1
    convert_bhv2mat_ML1(bhvdir,bhvfiles,pl2dir,lfpfile,...
        outdir,nicefile,psych_models);
elseif strcmp(filetype,'.bhv2') % ML2
    convert_bhv2mat_ML2(bhvdir,bhvfiles,pl2dir,lfpfile,...
        outdir,nicefile,psych_models);
end

close all

%% spit out number of trials in this session
load(strjoin({outdir,nicefile},'/'));

disp(['# free trials: ',num2str(sum(~isnan(trialinfo.lever) & trialinfo.lever~=0 & trialinfo.trialtype==2))])
disp(['# forced trials: ',num2str(sum(~isnan(trialinfo.lever) & trialinfo.lever~=0 & trialinfo.trialtype==1))])






