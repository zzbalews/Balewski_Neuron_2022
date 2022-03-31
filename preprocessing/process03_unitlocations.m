%% find sterotax coords for each unit

%% setup -- edit variables for session, then run!

addpath functions/neural/

addpath functions/'Matlab Offline Files SDK'/

dataHDDdir = '/media/StorageHDD/WallisLabData';

channelmap = '~/Documents/WallisLab/recording/Chap/Chap00_decoding_OFC+/Chap_info/channel_source.csv';
transform = 'Chap_gridcoord_transforms.mat';
% channelmap = '~/Documents/WallisLab/recording/George/George00_decoding_OFC+/George_info/channel_source.csv';
% transform = 'George_gridcoord_transforms.mat';

pl2_processed = 'pl2_processed';

outinfo = 'info';

%edit these for specific session
subject = 'Chap';
% experiment = 'decoding_CN+OFC';
% subject = 'George';
% experiment = 'decoding_ACC+OFC';
%%

for i = 1:13
    
    % chap files
    switch i
        case 1
            session_date = '11-29-2017';
            session_num = 'Chap00_rec11';
            experiment = 'decoding_ACC+OFC';
        case 2
            session_date = '12-11-2017';
            session_num = 'Chap00_rec16';
            experiment = 'decoding_ACC+OFC';
        case 3
            session_date = '11-28-2018';
            session_num = 'Chap00_rec60';
            experiment = 'decoding_ACC+OFC';
        case 4
            session_date = '11-30-2018';
            session_num = 'Chap00_rec61';
            experiment = 'decoding_ACC+OFC';
        case 5
            session_date = '12-10-2018';
            session_num = 'Chap00_rec62';
            experiment = 'decoding_ACC+OFC';
        case 6
            session_date = '12-14-2018';
            session_num = 'Chap00_rec63';
            experiment = 'decoding_ACC+OFC';
        case 7
            session_date = '01-02-2019';
            session_num = 'Chap00_rec65';
            experiment = 'decoding_ACC+OFC';
        case 8
            session_date = '01-08-2019';
            session_num = 'Chap00_rec66';
            experiment = 'decoding_ACC+OFC';
        case 9
            session_date = '04-17-2018';
            session_num = 'Chap00_rec33';
            experiment = 'decoding_DLPFC+OFC';
        case 10
            session_date = '05-04-2018';
            session_num = 'Chap00_rec36';
            experiment = 'decoding_DLPFC+OFC';
        case 11
            session_date = '06-13-2018';
            session_num = 'Chap00_rec41';
            experiment = 'decoding_DLPFC+OFC';
        case 12
            session_date = '06-15-2018';
            session_num = 'Chap00_rec42';
            experiment = 'decoding_DLPFC+OFC';
        case 13
            session_date = '11-22-2017';
            session_num = 'Chap00_rec09';
            experiment = 'decoding_DLPFC+OFC';
    end
    
    
    %     % george files
    %     switch i
    %         case 1
    %             session_date = '04-01-2021';
    %             session_num = 'George00_rec28';
    %         case 2
    %             session_date = '04-05-2021';
    %             session_num = 'George00_rec29';
    %         case 3
    %             session_date = '04-09-2021';
    %             session_num = 'George00_rec31';
    %         case 4
    %             session_date = '04-11-2021';
    %             session_num = 'George00_rec32';
    %         case 5
    %             session_date = '04-14-2021';
    %             session_num = 'George00_rec33';
    %     end
    
    %% get filenames
    
    session_name = strjoin({session_num,strrep(session_date,'-','')},'_');
    
    % spk data
    pl2_dir = strjoin({dataHDDdir,subject,experiment,pl2_processed,session_name},'/');
    temp = dir([pl2_dir,'/*-pics-SPKfr.mat']);
    spkfile = temp(1).name;
    
    % full out file
    out_dir = strjoin({dataHDDdir,subject,experiment,outinfo,session_name},'/');
    outfile = fullfile(out_dir,[session_name,'_unit_MLAPDVloc.mat']);
    
    
    %% load recording info & probe mapping
    recinfo = get_recording_info(out_dir,session_name);
    [~,probeID] = get_channel_regions(channelmap,session_date,session_num);
    
    %% get (ML,AP,DV) for each contact
    contactID = get_contact_location(subject,...
        probeID,recinfo,transform);
    
    %% assign (ML,AP,DV) loc to each unit
    save_unit_location(contactID,fullfile(pl2_dir,spkfile),outfile);
    
    
end
