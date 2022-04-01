function display_unit_facts(glmdata)


region_names = unique(glmdata.region);
subject_names = unique(glmdata.subject);

nreg = length(region_names);
nsubj = length(subject_names);


for R = 1:nreg
    for S = 1:nsubj
        
        % restrict units
        region = region_names{R};
        subject = subject_names{S};
                
        idx = strcmp(glmdata.region,region_names{R}) & ...
            strcmp(glmdata.subject,subject_names{S});
        
        disp([subject,' @ ',region])
        disp(['... # units = ',num2str(sum(idx))])
        
        % get total unique LFP channels (only with units)
        sessions = unique(glmdata.session(idx));
        track_ch = nan(length(sessions),1);
        for ses = 1:length(sessions)
           idx_ses = idx & strcmp(glmdata.session,sessions{ses});
           if sum(idx_ses) ~= 0
           ch = cellfun(@(x) str2num(x(9:11)),glmdata.unit_names(idx_ses));
           track_ch(ses) = length(unique(ch));
           end
        end
        disp(['... # channels = ',num2str(sum(track_ch))])
    end
end
end

