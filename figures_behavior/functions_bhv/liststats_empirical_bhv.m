function liststats_empirical_bhv(bhvdata)


subject_names = unique(bhvdata.subject);

temp = diff(bhvdata.subjval_expval,[],2);
bhvdata.high_opt = nan(size(temp));
bhvdata.high_opt(temp>0) = 1;
bhvdata.high_opt(temp<0) = -1;

for s = 1:length(subject_names)
    
    subject = subject_names{s};
    disp(subject)
    
    % # sessions
    nses = length(unique(bhvdata.session(strcmp(bhvdata.subject,subject))));
    disp(['... # sessions = ',num2str(nses)])
    
    % # trials
    tr_free = strcmp(bhvdata.subject,subject) & ...
        bhvdata.trialtype==2 & ...
        bhvdata.lever~=0;
    disp(['... avg # free = ',num2str(sum(tr_free)/nses)])
    tr_forced = strcmp(bhvdata.subject,subject) & ...
        bhvdata.trialtype==1 & ...
        bhvdata.lever~=0;
    disp(['... avg # forced = ',num2str(sum(tr_forced)/nses)])
    disp(['... avg # all trials = ',num2str(sum(tr_free+tr_forced)/nses)])
    
    % % chose high>low value?
    chose_high = bhvdata.lever(tr_free) == bhvdata.high_opt(tr_free);
    disp(['... % chose high>low free = ',num2str(100*mean(chose_high))])
    
    % avg RT
    disp(['... med RT free (ms) = ',num2str(median(bhvdata.rt(tr_free)))])
    disp(['... stdev RT free (ms) = ',num2str(std(bhvdata.rt(tr_free)))])
    disp(['... med RT forced (ms) = ',num2str(median(bhvdata.rt(tr_forced)))])
    disp(['... stdev RT forced (ms) = ',num2str(std(bhvdata.rt(tr_forced)))])
    
    % avg sacc 1 time
    disp(['... med sacc1 time free (ms) = ',num2str(nanmedian(bhvdata.MLsacc(tr_free,1)))])
    disp(['... stdev sacc1 time free (ms) = ',num2str(nanstd(bhvdata.MLsacc(tr_free,1))/sqrt(sum(tr_free)))])
    disp(['... med sacc1 time forced (ms) = ',num2str(nanmedian(bhvdata.MLsacc(tr_forced,1)))])
    disp(['... stdev sacc1 time forced (ms) = ',num2str(nanstd(bhvdata.MLsacc(tr_forced,1)))])
    
end

    
end

