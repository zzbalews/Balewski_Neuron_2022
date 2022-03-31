function make_lfp_cut(bhvdir,bhvfile,pl2dir,lfpfile,pl2processed_dir,...
    pics_range,choice_range,sacc1_range,sacc2_range,session_name,channels_with_spk)


%% determine file names, skip if output exist

% any LFP filtered files exist?
do_LFPs = 1;
outfile_prefix = strjoin({pl2processed_dir,[session_name,'-LFP-']},'/');
temp = dir([outfile_prefix,'*']);
if ~isempty(temp)
    disp(['LFP processing: ',num2str(length(temp)),' filtered files exist; STOP'])
    do_LFPs=0;
end

% thetaphase dir exists?
thetadir = strjoin({pl2processed_dir,'thetaphase'},'/');
if ~exist(thetadir)
    mkdir(thetadir)
end
theta_prefix = strjoin({thetadir,[session_name,'-LFP-02theta-PHASE-']},'/');


%% process LFP data

if do_LFPs
    %% extract LFP channel --> notch filter 60Hz --> bandpass --> cut into trials --> save output
    fname = strjoin({pl2dir,lfpfile},'/');
    
    % get channel names
    [~,names] = plx_adchan_names(fname);
    names = cellstr(names);
    names = names(cellfun(@(x) ~isempty(strfind(x,'FP')),names)); % only keep FP channels
    Nchannels = length(names);
    
    % get time properties
    [adfreq,nt] = plx_ad(fname,0);
    
    % build filters
    [notch_filts,band_filts,bandnames] = make_filters(adfreq,'.');
    Nnotch = length(notch_filts);
    Nband = length(band_filts);
    
    
    % adjust time to account for filter lag
    LFPtime = 1:nt; %ms
    lag = 1+(Nnotch+1)*adfreq/2;
    
    % get trial time cuts
    bhv = strjoin({bhvdir,bhvfile},'/');
    [pics_entries, choice_entries,sacc1_entries, sacc2_entries] = ...
        extract_trial_cuts(pics_range,choice_range,sacc1_range,sacc2_range,bhv,(nt-lag+2));
    
    tic
    % apply filters & cut trials (all channels)
    parfor ch = 1:Nchannels
        
        if ~ismember(str2num(names{ch}(3:5)),channels_with_spk)
            disp(['skipping   ...', names{ch}])
            continue
        end
        
        disp(['working on ...', names{ch}])
        
        % load LFP and apply notch filters
        working_signal = load_notched_LFP(fname,names{ch},notch_filts,'.',session_name);
        
        % for each band pass filter
        for b = 1:Nband
            
            % apply band filter and get magnitude/phase --> boxcar smooth --> zscore
            [mag, phs] = after_bandpass_LFP(working_signal, band_filts{b}.Coefficients,lag);
            
            % cut by trial
            pics_mag = cut_trials(mag,pics_entries);
            choice_mag = cut_trials(mag,choice_entries);
            sacc1_mag = cut_trials(mag,sacc1_entries);
            sacc2_mag = cut_trials(mag,sacc2_entries);
            
            % save
            outfile = [outfile_prefix,bandnames{b},'-',names{ch},'.mat'];
            save_mag(outfile,pics_mag,choice_mag,sacc1_mag,sacc2_mag,...
                pics_range,choice_range,sacc1_range,sacc2_range);
            
            
            % for theta only, save phase
            if strcmp(bandnames{b},'02theta')
                % cut by trial
                pics_phase = cut_trials(phs,pics_entries);
                choice_phase = cut_trials(phs,choice_entries);
                sacc1_phase = cut_trials(phs,sacc1_entries);
                sacc2_phase = cut_trials(phs,sacc2_entries);
                
                
                % save
                outfile = [theta_prefix,names{ch},'.mat'];
                save_phase(outfile,pics_phase,choice_phase,sacc1_phase,sacc2_phase,...
                    pics_range,choice_range,sacc1_range,sacc2_range);
            end
        end
    end
    toc
    
end
end