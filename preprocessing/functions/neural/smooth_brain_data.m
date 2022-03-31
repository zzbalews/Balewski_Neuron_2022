function smooth_brain_data(pl2processed_dir, pl2cut_dir, channelID, smth_sets, session_name)

%% init data files or skip if already saved

% check for smoothing insturctions
if isempty(fieldnames(smth_sets))
    error('missing smoothing instructions! use add_smoothing')
end

% init missing files
[fnames,needed_regions,needed_smooths,needed_trialperiods] = identify_missing_files(smth_sets,session_name,pl2cut_dir);

if isempty(needed_regions) || isempty(needed_smooths)
    disp('all smoothing is done!')
    return
end

%% load SPKfr data
disp('working on smoothing ... SPKfr')


for T = 1:length(needed_trialperiods)
    
    % load units, synced to trial period
    [spkfr_allunits, spkfr_range, unit_names] = loadspk(pl2processed_dir,'SPKfr',needed_trialperiods{T});
    
    if T==1 % just once, get useful info for all smoothings
        
        [ntr,ntimes,~] = size(spkfr_allunits);
        
        % assign units to regions
        [~,region_units] = get_region_index(unit_names,channelID,needed_regions,'spk');
        
        % find midpoints for smoothing
        [timestamp, midpoints] = get_smth_points(smth_sets,needed_smooths,spkfr_range);
        
    end
    
    % do smoothing
    for i = 1:length(needed_smooths)
        
        % apply fastsmooth
        iter = needed_smooths(i);
        spkfr_allunits_smooth = smooth_array(spkfr_allunits, timestamp, midpoints{i}, smth_sets(iter).w);
        t_mids = midpoints{i};
        
        % split by region
        for R = 1:length(needed_regions)
            
            % save if file is missing
            if ~isempty(fnames(iter).(needed_trialperiods{T}).(needed_regions{R}))
                SPKfr_units = spkfr_allunits_smooth(:,:,region_units{R});
                SPKfr_names = unit_names(region_units{R});
                
                save(strjoin({pl2cut_dir,fnames(iter).(needed_trialperiods{T}).(needed_regions{R})},'/'),...
                    'SPKfr_units','SPKfr_names','t_mids','-append');
                
                clear SPKfr_units 
            end
        end
        
        clear SPKfr_allunits_smooth
    end
    
    clear spkfr_allunits
end


%% load SPKfrnorm data
disp('working on smoothing ... SPKfrnorm')


for T = 1:length(needed_trialperiods)
    
    % load units, synced to trial period
    [spkfrnorm_allunits, spkfrnorm_range, unit_names] = loadspk(pl2processed_dir,'SPKfrnorm',needed_trialperiods{T});
    if spkfr_range~=spkfrnorm_range
        error('SPKfr and SPKfrnorm: different window cuts')
    end
    
    % do smoothing
    for i = 1:length(needed_smooths)
        
        % apply fastsmooth
        iter = needed_smooths(i);
        spkfr_allunits_smooth = smooth_array(spkfrnorm_allunits, timestamp, midpoints{i}, smth_sets(iter).w);
        t_mids = midpoints{i};
        
        % split by region
        for R = 1:length(needed_regions)
            
            % save if file is missing
            if ~isempty(fnames(iter).(needed_trialperiods{T}).(needed_regions{R}))
                SPKfrnorm_units = spkfr_allunits_smooth(:,:,region_units{R});
                SPKfrnorm_names = unit_names(region_units{R});
                
                save(strjoin({pl2cut_dir,fnames(iter).(needed_trialperiods{T}).(needed_regions{R})},'/'),...
                    'SPKfrnorm_units','SPKfrnorm_names','-append');
                
                clear SPKfrnorm_units
            end
        end
        
        clear spkfr_allunits_smooth
    end
    
    clear spkfrnorm_allunits
end

%% display # of units
for R = 1:length(needed_regions)
    disp(['# units in ',needed_regions{R},': ',num2str(sum(region_units{R}))]);
end

%% load LFP files
freqbands = {'01delta','02theta','03alpha','04beta','05gamma','06highgamma'};

for freq = 1:length(freqbands)
    
    disp(['working on smoothing ... ',freqbands{freq}])
    
    % get file names
    temp = dir(strjoin({pl2processed_dir,['*LFP-',freqbands{freq},'*']},'/'));
    if isempty(temp)
        error('check LFP files: first band missing')
    end
    freqfiles = {temp.name};
    
    
    if freq==1 % for first freq, get channel regions
        [channelnames,region_channels] = get_region_index(freqfiles,channelID,needed_regions,'lfp');
    else % for remaining, ensure same number of channels
        if length(freqfiles)~=length(channelnames)
            error('check LFP files: different number of channels in freq bands')
        end
    end
    
    % load, smooth, and save
    for R = 1:length(needed_regions) % for each region
        
        for T = 1:length(needed_trialperiods) % for each trial period
            
            % combine all channels in region
            idx = find(region_channels{R});
            lfp_allchannels = nan(ntr,ntimes,length(idx));
            parfor f = 1:length(idx)
                [lfp_range,lfp_mag] = loadlfp(strjoin({pl2processed_dir,freqfiles{idx(f)}},'/'),needed_trialperiods{T});
                lfp_allchannels(:,:,f) = lfp_mag;
                
                if lfp_range~=spkfr_range
                    error('SPKfr and LFP do not have same cut window')
                end
            end
            
            % do smoothing
            for i = 1:length(needed_smooths)
                iter = needed_smooths(i);
                
                
                if freq==1 % for first band, save channel names
                    LFP_names = channelnames(idx);
                    save(strjoin({pl2cut_dir,fnames(iter).(needed_trialperiods{T}).(needed_regions{R})},'/'),...
                        'LFP_names','-append');
                end
                
                % apply fastsmooth
                eval(['LFP',freqbands{freq},'= smooth_array(lfp_allchannels,timestamp,midpoints{i}, smth_sets(iter).w);'])
                
                % save
                save(strjoin({pl2cut_dir,fnames(iter).(needed_trialperiods{T}).(needed_regions{R})},'/'),...
                    ['LFP',freqbands{freq}],'-append');
                
            end
        end
    end
end
end