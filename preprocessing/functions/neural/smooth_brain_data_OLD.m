function ses_info = smooth_brain_data(ses_info, regions, smoothtype, w, s, slices, dir_pl2_smooth)
% for breain region = reg
% if smoothtype='slice', avg spk and lfp at t=slices (1x2 array for pics,
% choice) using window=w
% if smoothtype='boxcar', smooth spk and lfp with boxcar width w and step s;
% cat into one matrix and save in dir_pl2_smooth

%% setup
% make filename(s) for regions
if strcmp(smoothtype,'slice') % one slice for free decoding
    outfile_prefix = [dir_pl2_smooth, ses_info.name, ...
        '-cutbraindata-slice',num2str(w),'ms-'];
else %if strcmp(smoothtype,'boxcar') % boxcar sliding avg
    outfile_prefix = [dir_pl2_smooth, ses_info.name, ...
        '-cutbraindata-w',num2str(w),'-s',num2str(s),'-'];
end

disp(['working on ',outfile_prefix,'*'])

% if exists, skip!
missing_regions = ones(size(regions));
for r = 1:length(regions)
    ses_info.(regions{r}).(smoothtype) = [outfile_prefix,regions{r},'-'];
    
    
    temp = dir([ses_info.(regions{r}).(smoothtype),'*']);
    if length(temp)>0
        disp([regions{r}, ' smoothed data saved!'])
        missing_regions(r)=0;
    end
end

if sum(missing_regions)==0
    return
end

brain_data = struct();

%% get time midpoints
load(ses_info.spkfile,'pics_range','choice_range');

[t_pics, t_mids_pics, t_choice, t_mids_choice] = get_midpoints(smoothtype, ...
    slices, pics_range, choice_range, w, s);

%% smooth spk data (split by trial period)
disp('...smoothing spk');

% load spk
load(ses_info.spkfile,'pics_allunits','choice_allunits','unit_names');

for r = find(missing_regions)
    % keep only region units
    chan = cellfun(@(x) str2num(x(9:11)), unit_names);
    keepunits = find(ismember(chan, ses_info.(regions{r}).channels));
    
    % save channel and code spikes as freq band = 0
    brain_data.(regions{r}).spk_channel = chan(keepunits);
    brain_data.(regions{r}).spk_freq = zeros(size(brain_data.(regions{r}).spk_channel));
    
    % smooth arrays
    data = smooth_array(pics_allunits(:,:,keepunits), ...
        t_pics, t_mids_pics, w);
    brain_data.(regions{r}).pics_allunits_smooth = data;
    [ntr,nmids_pics,~] = size(pics_allunits);
    
    
%     data = smooth_array(choice_allunits(:,:,keepunits), ...
%         t_choice, t_mids_choice, w);
%     brain_data.(regions{r}).choice_allunits_smooth = data;
    
end

%% smooth lfp data (split by trial period)
disp('...smoothing lfp');

for r = find(missing_regions)
    % keep only region channels
    chan = cellfun(@(x) str2num(x(end-6:end-4)), ses_info.lfpfile);
    freq = cellfun(@(x) str2num(x(strfind(x,'-LFP-')+(5:6))), ses_info.lfpfile);
    keeplfp = find(ismember(chan, ses_info.(regions{r}).channels));
    
    % save channel and freq band
    brain_data.(regions{r}).lfp_channel = chan(keeplfp);
    brain_data.(regions{r}).lfp_freq = freq(keeplfp);
    
    % load all lfp
    [pics_allmags, choice_allmags] = load_lfpfiles(ses_info.lfpfile(keeplfp),ntr,nmids_pics);
    
    % smooth arrays
    data = smooth_array(pics_allmags,...
        t_pics, t_mids_pics, w);
    brain_data.(regions{r}).pics_allmags_smooth = data;
    
%     choice_allmags_smooth = smooth_array(choice_allmags,...
%         t_choice, t_mids_choice, w);
%     brain_data.(regions{r}).choice_allmags_smooth = data;
    
end
%% combine spk and lfp; save
disp('...combine and save');

for r = find(missing_regions)
    pics = cat(3, brain_data.(regions{r}).pics_allunits_smooth, ...
        brain_data.(regions{r}).pics_allmags_smooth);
    
%     choice = cat(3, brain_data.(regions{r}).choice_allunits_smooth, ...
%         brain_data.(regions{r}).choice_allmags_smooth);
    
    channels = [brain_data.(regions{r}).spk_channel; ...
        brain_data.(regions{r}).lfp_channel];
    freqbands = [brain_data.(regions{r}).spk_freq; ...
        brain_data.(regions{r}).lfp_freq];
    
    save([ses_info.(regions{r}).(smoothtype),'pics.mat'], ...
        'pics','channels','freqbands','t_mids_pics','-v7.3');
%     save([ses_info.(regions{r}).(smoothtype),'choice.mat'], ...
%         'choice','channels','freqbands','t_mids_choice','-v7.3');

end

end

