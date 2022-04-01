function out_states = get_states(postprob_delta, nlevels, nconsec, thresh)
% threshold = 2x baseline or 200% (1.95x or 195% if nlevels = 2); get
% states from postprob_delta; state must be >=nconsec bins long


ntr = size(postprob_delta,1);

out_states = nan([size(postprob_delta),length(thresh)]);

for i = 1:length(thresh)
    above_thresh = postprob_delta >=thresh(i);
    
    for tr = 1:ntr
        for lev = 1:nlevels
            out_states(tr,:,:,lev,i) = give_consec_seg(above_thresh(tr,:,:,lev), nconsec);
        end
    end
end

end

