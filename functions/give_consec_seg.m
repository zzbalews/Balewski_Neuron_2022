function [ track_sig_consec,output ] = give_consec_seg( track_sig, N )

track_sig_consec = zeros(size(track_sig));

h_idx = find(track_sig);
h_diff = [10 diff(h_idx)];

seg_starts = find(h_diff>1);
if h_diff(end)==1
    seg_starts(end+1) = length(h_diff)+1;
else
    seg_starts(end+1) = seg_starts(end)+1;
end

seg_lens = diff(seg_starts);

keep = find(seg_lens(1:end)>=N);


for k = 1:length(keep)
    idx = h_idx(seg_starts(keep(k))) + (0:seg_lens(keep(k))-1) ;
    track_sig_consec(idx)=1;
end

output = seg_lens(keep);



end

