function [channelID, probeID] = get_channel_regions(channelmap,session_date,session_num)

% load all sessions
locations = readtable(channelmap,'FileType','delimitedtext','Delimiter',',');

% % ignore L/R hemi divion?
% locations(:,4) = cellfun(@(x) x(2:end),locations{:,4},'UniformOutput',0);

% convert date
nicedate = strjoin(strsplit(session_date,'-'),'/');

% pull channels for this session
locations_ses = locations(strcmp(locations.SessionName,session_num),:);

% get regions
regions = unique(locations_ses.Region);

% save channels for each region
channelID = struct();
for r = 1:length(regions)
    idx = strcmp(locations_ses.Region,regions{r});
   channelID.(regions{r}) = locations_ses(idx,'ChannelNumber');
end

probeID = locations_ses(:,3:5);

end

