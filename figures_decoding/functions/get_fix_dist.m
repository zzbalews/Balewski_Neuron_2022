function [dirdata_pics,dirdata_choice] = get_fix_dist(bhvdata, dirdata_pics, dirdata_choice)
%%
subject_names = unique(bhvdata.subject);
nsubj = length(subject_names);

t_mids = dirdata_pics.t_mids;
nmids = length(t_mids);


t_fix = t_mids>=-500 & t_mids<0;
t_pics = t_mids>=0 & t_mids<500;
t_choice = t_mids>=-500 & t_mids<0;

if isempty(dirdata_choice)
    skipchoice = 1;
else
    skipchoice = 0;
end

dirdata_pics.postprob_ch_disc = nan(size(dirdata_pics.postprob_ch));
dirdata_pics.postprob_ch_norm = nan(size(dirdata_pics.postprob_ch));

if skipchoice
    dirdata_choice = [];
else
    dirdata_choice.postprob_ch_disc = nan(size(dirdata_choice.postprob_ch));
    dirdata_choice.postprob_ch_norm = nan(size(dirdata_choice.postprob_ch));
end

for S = 1:nsubj
    
    subject = subject_names{S};
    idx_subj = strcmp(bhvdata.subject,subject);
    sessions = unique(bhvdata.session(idx_subj,:));
    
    nses = length(sessions);
    track_quant = nan(nses,2);
    track_m = nan(nses,1);
    track_sd = nan(nses,1);
    
    for ses = 1:nses
        
        % get trials for session
        idx_ses = strcmp(bhvdata.session,sessions{ses});
        ntr = sum(idx_ses);
        
        postprob_ses = dirdata_pics.postprob_ch(idx_ses,:);
        
        % build norm dist from fix period
        fixdist = reshape(postprob_ses(:,t_fix),1,[]);
        track_quant(ses,:) = quantile(fixdist,[.01 .99]);
        %         track_quant(ses,:) = quantile(fixdist,[.05 .95]);
        
        track_m(ses) = nanmean(fixdist);
        track_sd(ses) = nanstd(fixdist);
        
        % norm pics by fix dist & get states
        disc_ses = zeros(ntr,nmids);
        ch = postprob_ses>=track_quant(ses,2);
        unch = postprob_ses<=track_quant(ses,1);
        
        for tr = 1:ntr
            ch(tr,:) = give_consec_seg(ch(tr,:),4);
            unch(tr,:) = give_consec_seg(unch(tr,:),4);
        end
        
        disc_ses(ch) = 1;
        disc_ses(unch) = -1;
        
        dirdata_pics.postprob_ch_disc(idx_ses,:) = disc_ses;
        dirdata_pics.postprob_ch_norm(idx_ses,:) = (postprob_ses-track_m(ses))./track_sd(ses);
        
        if ~skipchoice
            % norm choice by fix dist & get states
            postprob_ses = dirdata_choice.postprob_ch(idx_ses,:);
            disc_ses = zeros(ntr,nmids);
            
            ch = postprob_ses>=track_quant(ses,2);
            unch = postprob_ses<=track_quant(ses,1);
            
            for tr = 1:ntr
                ch(tr,:) = give_consec_seg(ch(tr,:),4);
                unch(tr,:) = give_consec_seg(unch(tr,:),4);
            end
            
            disc_ses(ch) = 1;
            disc_ses(unch) = -1;
            
            dirdata_choice.postprob_ch_disc(idx_ses,:) = disc_ses;
            dirdata_choice.postprob_ch_norm(idx_ses,:) = (postprob_ses-track_m(ses))./track_sd(ses);
        end
        
    end
    temp = mean(mean(abs(track_quant-.5),2));
    disp([subject,': 99% cutoff=',num2str(0.5+round(temp,3))])
end


end
