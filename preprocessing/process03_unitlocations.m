%% find sterotax coords for each unit
% all preprocessing scripts run for this recording session: George rec18 02-11-2021
% swap out session details and file paths to run for other sessions


%% setup -- edit variables for session, then run!

addpath functions/neural/
addpath functions/'Matlab Offline Files SDK'/

% get location for each channel in recordings
channelmap = '../info/George_channel_source.csv'; % recording location for each channel
transform = '../info/George_gridcoord_transforms.mat';

data_dir = '../data_raw'; % where to find raw neural pl2 data
pl2_temp = '../data_temp'; % where to find mid processing neural data
out_dir = '../data_processed'; % where to save unit locations

% edit these for specific session
subject = 'George';
experiment = 'decoding_CN+OFC';
session_date = '02-11-2021';
session_num = 'George00_rec18';

%% get filenames
session_name = strjoin({session_num,strrep(session_date,'-','')},'_');

% spk data
spkfile = 'George00_rec18_02112021-pics-SPKfr.mat';

% full out file
outfile = fullfile(out_dir,[session_name,'_unit_MLAPDVloc.mat']);

% recording details file
recfile = fullfile(data_dir,'George00_rec18_02112021_elec_info_clean.csv');

%% load recording info & probe mapping
recinfo = get_recording_info(recfile);
[~,probeID] = get_channel_regions(channelmap,session_date,session_num);

%% get (ML,AP,DV) for each contact
contactID = get_contact_location(subject,...
    probeID,recinfo,transform);

%% assign (ML,AP,DV) loc to each unit
save_unit_location(contactID,fullfile(pl2_temp,spkfile),outfile);


