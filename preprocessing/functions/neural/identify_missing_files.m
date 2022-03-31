function [fnames,needed_regions,needed_smooths,needed_trialperiods] = identify_missing_files(smth_sets,session_name,pl2cut_dir)


% which trial periods and regions are asked for?
needed_trialperiods = unique(cat(2,smth_sets.trialperiods));
possible_regions = unique(cat(2,smth_sets.regions));

% which regions and smoothings options are still needed?
needed_regions = {};
needed_smooths = [];

% output file names, if still needed
fnames = struct();

for R = 1:length(possible_regions) % for all possible regions
    
    for T = 1:length(needed_trialperiods) % for all asked trial periods
        
        for iter = 1:length(smth_sets) % for all asked smoothigns
            
            % build file name
            thisfilename = [session_name,...
                smth_sets(iter).label,'_',needed_trialperiods{T},'_',possible_regions{R},'.mat'];
            
            % check if file exists
            if exist(strjoin({pl2cut_dir,thisfilename},'/')) % yes! 
                % don't do anything
                fnames(iter).(needed_trialperiods{T}).(possible_regions{R}) = []; 
                
            else % missing!
                % make file
                fnames(iter).(needed_trialperiods{T}).(possible_regions{R}) = thisfilename;
                
                meta = struct();
                meta.window = smth_sets(iter).w;
                meta.step = smth_sets(iter).s;
                meta.region = possible_regions{R};
                meta.trialperiod = needed_trialperiods{T};
                save(strjoin({pl2cut_dir,thisfilename},'/'),'meta','thisfilename','-v7.3')
                
                % add to to-do list
                needed_regions{end+1} = possible_regions{R};
                needed_smooths(end+1) = iter;
            end
            
        end
    end
end

needed_regions = unique(needed_regions);
needed_smooths = unique(needed_smooths);


end

