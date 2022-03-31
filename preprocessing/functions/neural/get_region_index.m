function [names,region_idx] = get_region_index(longnames,channelID,missing_regions,signal)

if strcmp(signal,'spk')
    % use unit_names
    names = longnames;
    names_numbersonly = cellfun(@(x) str2num(x(9:11)),names);
    
elseif strcmp(signal,'lfp')
% get shorter channel names from files
names = cellfun(@(x) x(end-8:end-4),longnames,'UniformOutput',false);
names_numbersonly = cellfun(@(x) str2num(x(3:5)),names);

end

% assign freqfile index to region
region_idx = cell(size(missing_regions));

for R = 1:length(missing_regions)
    temp = table2array(channelID.(missing_regions{R}));
    region_idx{R} =  ismember(names_numbersonly,temp);
end

end

