function checkout_decoder_valdircorr(bhvdata,decdata1,decdata2,keep_tr,which_level)
%%
% useful vars

subject_names = unique(bhvdata.subject);
nsubj = length(subject_names);

bounds = [-850,1200];

t_mids = decdata1.t_mids;
t_subset = find(t_mids>=bounds(1) & t_mids<bounds(2));
T = t_mids(t_subset);

niters = 1000;

pairs = [1,2; 1,3; 1,4;2,3;2,4;3,4];

if strcmp(which_level,'chosen+dir')
clr = [0 .7 0];
ppd_dir = decdata1.postprob_ch;
ppd_val = decdata2.ch_ppd;
show_avg = 1;  
end


for s = 1:nsubj
    subject = subject_names{s};
    idx = find(keep_tr & ...
        strcmp(bhvdata.subject, subject));
    
    ntr = length(idx);
    
    if ntr==0
        continue
    end
    
    mini_PPD_dir = ppd_dir(idx,t_subset)-.5;
    mini_PPD_val = ppd_val(idx,t_subset);
    
    rescale_dir = 1/max(mean(mini_PPD_dir));
    rescale_val = 1/max(mean(mini_PPD_val));
    
    mini_PPD_dir = mini_PPD_dir * rescale_dir;
    mini_PPD_val = mini_PPD_val * rescale_val;
    
    figure;
    if show_avg
    subplot(2,1,1); hold on
    plot(bounds,[0 0],':','Color',[.5 .5 .5])
    plot([0 0],[-.1 1.1],':','Color',[.5 .5 .5])
    plot([-750 -750],[-.1 1.1],':','Color',[.5 .5 .5])
    
    a=shadedErrorBar(T,mean(mini_PPD_dir),...
        std(mini_PPD_dir)./sqrt(ntr),...
        'lineprops',{'Color',[.7 0 .7]});
    b=shadedErrorBar(T,mean(mini_PPD_val),...
        std(mini_PPD_val)./sqrt(ntr),...
        'lineprops',{'Color',[0 .7 0]});
    ylim([-.1 1.1])
    xlim(bounds)
    xlabel('time from pics (ms)')
    ylabel([which_level,'normalized decoding'])
    legend([a.mainLine,b.mainLine],{'CN direction','CN value'},'box','off')
    title(subject)
    end

    
    
    subplot(2,1,2); hold on
    plot(bounds,[0 0],':','Color',[.5 .5 .5])
    plot([0 0],[-.15 .3],':','Color',[.5 .5 .5])
    plot([-750 -750],[-.15 .3],':','Color',[.5 .5 .5])
    
    ylim([-.15 .3])
    xlim(bounds)
    
    
    % show shuffled correlation distribution
    shufCorr_ch = nan(niters,length(T));
    
    for i = 1:niters
        if mod(i,100)==0
            disp(i)
        end
        row_shuf = nan(ntr,1);
        for p = 1:size(pairs,1)
            idx_pair = find(min(bhvdata.valbin_expval(idx,:),[],2)==pairs(p,1) & ...
                max(bhvdata.valbin_expval(idx,:),[],2)==pairs(p,2));
            row_shuf(idx_pair) = idx_pair(randperm(length(idx_pair)));
        end
        
%         row_shuf = randperm(ntr);
        OFC_PPD_ch_shuf = mini_PPD_dir(row_shuf,:);
        temp = corr(OFC_PPD_ch_shuf, mini_PPD_val);
        shufCorr_ch(i,:) = temp(eye(length(T))==1);
        
    end
    
     shuf_ch_low = quantile(shufCorr_ch,.01,1);
    shuf_ch_med = quantile(shufCorr_ch,.5,1);
    shuf_ch_high = quantile(shufCorr_ch,.99,1);
    
    b=shadedErrorBar(T,shuf_ch_med,...
        [shuf_ch_med-shuf_ch_low; shuf_ch_high-shuf_ch_med],...
        'lineprops',{'Color',[.6 .6 .6]});
    
    % show real correlations
    temp = corr(mini_PPD_dir,mini_PPD_val);
    rawCorr_ch = temp(eye(length(T))==1);
    a=plot(T,rawCorr_ch,'Color',[clr .5],'LineWidth',1);
    ifsig = double(rawCorr_ch>=shuf_ch_high');
    ifsig(ifsig==0) = NaN;
    plot(T,rawCorr_ch.*ifsig,'Color',clr,'LineWidth',2)
    
    
    legend([a,b.mainLine],{'real',['shuffled x',num2str(niters)]},'box','off')
    xlabel('time from pics (ms)')
    ylabel(['correlation OFC:CN (',which_level,' PPD)'])
    
    
    
    
end


end


