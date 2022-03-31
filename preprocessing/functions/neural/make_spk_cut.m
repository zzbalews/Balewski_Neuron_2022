function [channels_with_spk] = make_spk_cut(bhvdir,bhvfile,pl2dir,spkfiles,pl2processed_dir,...
    pics_range,choice_range,sacc1_range,sacc2_range,session_name,do_zscore,which_sync,restricted_channels)


%% determine file names, skip if output exist

if do_zscore % normalize firing rate
    
    outfile = strjoin({pl2processed_dir,[session_name,'-',which_sync,'-SPKfrnorm.mat']},'/');
    
else % raw firing rates
    
    outfile = strjoin({pl2processed_dir,[session_name,'-',which_sync,'-SPKfr.mat']},'/');
    
end

temp = {'fr','frnorm'};
do_cut = 1;
if exist(outfile)
    disp('SPK preprocessing: STOP')
    disp(['    (cut ',which_sync,' ',temp{do_zscore+1},' exists)'])
    do_cut = 0;
else
    disp('SPK preprocessing: continue')
    disp(['    (missing ',which_sync,' ',temp{do_zscore+1},')'])
    
end



%% process spike data

if do_cut
    %% extract spike times --> spike trains --> smooth & norm
    
    % get number spk files
    nspk = length(spkfiles);
    
    % get time cuts
    [pics_entries, choice_entries, sacc1_entries, sacc2_entries] = extract_trial_cuts(...
        pics_range,choice_range,sacc1_range,sacc2_range,strjoin({bhvdir,bhvfile},'/'),NaN);
    
    % max time required
    maxT = max([reshape(pics_entries,1,[]),reshape(choice_entries,1,[])]);
    
    % save unique unit names
    unit_names = {};
    
    % save unit firing rates
    raw_firing_rate = [];
    
    % save cut data
    if strcmp(which_sync,'pics')
        pics_allunits= nan([size(pics_entries),0]);
    elseif strcmp(which_sync,'choice')
        choice_allunits = nan([size(choice_entries),0]);
    elseif strcmp(which_sync,'sacc1')
        sacc1_allunits = nan([size(sacc1_entries),0]);
    elseif strcmp(which_sync,'sacc2')
        sacc2_allunits = nan([size(sacc1_entries),0]);
    end
    
    
    % for all unit files associated with session
    total_channels = 0;
    for f = 1:nspk
        % get pl2 filename
        spkfname = strjoin({pl2dir,spkfiles{f}},'/');
        
        % get channel names
        [Nchannels,names] = plx_chan_names(spkfname);
        names = cellstr(names);
        
        % only look at restricted channels
        pull_channels = find(ismember(cellfun(@(x) str2num(x(9:11)), names), restricted_channels));
        if isempty(pull_channels)
            continue
        end
        
        % load timestamps of all waveforms of labeled units
        for ch = pull_channels'
            
            for u = 1:6 % up to 6 units per channel (never gets this high)
                
                % get time stamps
                [~,ts] = plx_ts(spkfname,ch,u);
                
                if ts == -1 % exceeded # units on this channel, continue
                    break
                else
                    % from spike times, get train
                    [spiketrain_smooth, unit_names, raw_fr] = get_unit_train(do_zscore,...
                        ts, unit_names, names{ch}, u, maxT);
                    
                    % nth unit
                    unitN = length(unit_names);
                    
                    % update raw firing rate
                    raw_firing_rate(unitN) = raw_fr; %Hz
                    
                    % cut data
                    if strcmp(which_sync,'pics')
                        pics_fr = cut_trials(spiketrain_smooth,pics_entries);
                        pics_allunits(:,:,unitN) = pics_fr;
                    elseif strcmp(which_sync,'choice')
                        choice_fr = cut_trials(spiketrain_smooth,choice_entries);
                        choice_allunits(:,:,unitN) = choice_fr;
                    elseif strcmp(which_sync,'sacc1')
                        sacc1_fr = cut_trials(spiketrain_smooth,sacc1_entries);
                        sacc1_allunits(:,:,unitN) = sacc1_fr;
                    elseif strcmp(which_sync,'sacc2')
                        sacc2_fr = cut_trials(spiketrain_smooth,sacc2_entries);
                        sacc2_allunits(:,:,unitN) = sacc2_fr;
                        
                    end
                    
                end
            end
        end
        
        total_channels = total_channels + length(unique(names));
    end
    
    
    %% save files
    if strcmp(which_sync,'pics')
        
        save(outfile, 'pics_allunits',...
            'pics_range','unit_names','raw_firing_rate','-v7.3')
        
    elseif strcmp(which_sync,'choice')
        save(outfile, 'choice_allunits',...
            'choice_range','unit_names','raw_firing_rate','-v7.3')
        
    elseif strcmp(which_sync,'sacc1')
        save(outfile, 'sacc1_allunits',...
            'sacc1_range','unit_names','raw_firing_rate','-v7.3')
        
    elseif strcmp(which_sync,'sacc2')
        save(outfile, 'sacc2_allunits',...
            'sacc2_range','unit_names','raw_firing_rate','-v7.3')
        
    end
    
    
    % get channels with units
    channels_with_spk = unique(cellfun(@(x) str2num(x(9:11)),unit_names));
    
    disp(['# unique unit channels: ',num2str(total_channels)]);
    
    
else % get channels with units
    load(outfile,'unit_names');
    channels_with_spk = unique(cellfun(@(x) str2num(x(9:11)),unit_names));
end
end