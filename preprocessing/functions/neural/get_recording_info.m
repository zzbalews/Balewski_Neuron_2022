function [recinfo] = get_recording_info(out_dir,session_name)

% load recording info file
temp = dir([out_dir,'/*elec_info_clean.csv']);
fname = fullfile(out_dir,temp(1).name);

recinfo = readtable(fname,'FileType','delimitedtext','Delimiter',',');


end

