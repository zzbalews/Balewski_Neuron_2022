function [recinfo] = get_recording_info(recfile)

% load recording info file

recinfo = readtable(recfile,'FileType','delimitedtext','Delimiter',',');


end

