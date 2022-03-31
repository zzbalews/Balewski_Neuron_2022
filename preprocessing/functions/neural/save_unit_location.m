function save_unit_location(contactID,spkfile,outfile)

% load spk unit names
load(spkfile,'unit_names','raw_firing_rate');

% extract channel from unit name
unit_channel = cellfun(@(x) str2num(x(9:11)), unit_names);
nunits = length(unit_channel);

% lookup (ML,AP,DV) for each unit
coords_ML_AP_DV = nan(nunits,3);
region = cell(nunits,1);

for u = 1:nunits
    idx = find(contactID{:,1}==unit_channel(u));
    
    coords_ML_AP_DV(u,:) = contactID{idx,{'ML','AP','DV'}};

    region{u} = contactID{idx,'Region'}{1};
end

% save!
save(outfile,'unit_names','raw_firing_rate','coords_ML_AP_DV','region');

end

