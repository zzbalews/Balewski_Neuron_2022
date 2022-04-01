function checkout_decoder_correlations(bhvdata,decdata1,decdata2,keep_tr,which_level)
%%
% useful vars

subject_names = unique(bhvdata.subject);
nsubj = length(subject_names);

bounds = [-950 950];

t_mids = decdata1.t_mids;
t_subset = find(t_mids>=bounds(1) & t_mids<bounds(2));
T = t_mids(t_subset);

niters = 10000;

pairs = [1,2; 1,3; 1,4;2,3;2,4;3,4];

if strcmp(which_level,'chosen')
clr = [1 0 0];
ppd_OFC = decdata1.ch_ppd;
ppd_CN = decdata2.ch_ppd;
show_avg = 1;
elseif strcmp(which_level,'unchosen')
    clr = [0 0 1];
   ppd_OFC = decdata1.unch_ppd;
ppd_CN = decdata2.unch_ppd;
   show_avg = 1;
 
elseif strcmp(which_level,'states')
    clr = [1 0 1];
     ppd_OFC = decdata1.ch_state - decdata1.unch_state;
ppd_CN = decdata2.ch_state - decdata2.unch_state;
  show_avg = 0;
elseif strcmp(which_level,'val+dir')
    clr = [0 .8 0];
    ppd_OFC = decdata1.ch_ppd;
    ppd_CN = decdata2.postprob_ch;
  
end
figure;

for s = 1:nsubj
    subject = subject_names{s};
    idx = find(keep_tr & ...
        strcmp(bhvdata.subject, subject));
    
    ntr = length(idx);
    
    if ntr==0
        continue
    end
    
    OFC_PPD_ch = ppd_OFC(idx,t_subset);
    CN_PPD_ch = ppd_CN(idx,t_subset);
    
    
%     figure;
%     if show_avg
%     subplot(2,1,1); hold on
%     plot(bounds,[0 0],':','Color',[.5 .5 .5])
%     plot([0 0],[-40 150],':','Color',[.5 .5 .5])
%     plot([-750 -750],[-40 150],':','Color',[.5 .5 .5])
%     
%     a=shadedErrorBar(T,mean(OFC_PPD_ch),...
%         std(OFC_PPD_ch)./sqrt(ntr),...
%         'lineprops',{'Color',[.7 0 .7]});
%     b=shadedErrorBar(T,mean(CN_PPD_ch),...
%         std(CN_PPD_ch)./sqrt(ntr),...
%         'lineprops',{'Color',[0 .7 0]});
%     ylim([-40 150])
%     xlim(bounds)
%     xlabel('time from pics (ms)')
%     ylabel([which_level,'PPD (value decoding)'])
%     legend([a.mainLine,b.mainLine],{'OFC','CN'},'box','off')
%     title(subject)
%     end
% 
%     subplot(2,1,2); 
    subplot(1,2,s);
hold on
    plot(bounds,[0 0],':','Color',[.5 .5 .5])
    plot([0 0],[-.15 .3],':','Color',[.5 .5 .5])
    plot([-750 -750],[-.15 .3],':','Color',[.5 .5 .5])
    
    ylim([-.15 .3])
    xlim(bounds)
    
    
    % show shuffled correlation distribution
    shufCorr_ch = nan(niters,length(T));
    
    for i = 1:niters
        if mod(i,100)==0
            disp([num2str(i),' / ',num2str(niters)])
        end
        row_shuf = nan(ntr,1);
        for p = 1:size(pairs,1)
            idx_pair = find(min(bhvdata.valbin_expval(idx,:),[],2)==pairs(p,1) & ...
                max(bhvdata.valbin_expval(idx,:),[],2)==pairs(p,2));
            row_shuf(idx_pair) = idx_pair(randperm(length(idx_pair)));
        end
        
%         row_shuf = randperm(ntr);
        OFC_PPD_ch_shuf = OFC_PPD_ch(row_shuf,:);
        temp = corr(OFC_PPD_ch_shuf, CN_PPD_ch);
        shufCorr_ch(i,:) = temp(eye(length(T))==1);
        
    end
    
     shuf_ch_low = quantile(shufCorr_ch,.01,1);
    shuf_ch_med = quantile(shufCorr_ch,.5,1);
    shuf_ch_high = quantile(shufCorr_ch,.99,1);
    
    b=shadedErrorBar(T,shuf_ch_med,...
        [shuf_ch_med-shuf_ch_low; shuf_ch_high-shuf_ch_med],...
        'lineprops',{'Color',[.6 .6 .6]});
    
    % show real correlations
    temp = corr(OFC_PPD_ch,CN_PPD_ch);
    rawCorr_ch = temp(eye(length(T))==1);
    a=plot(T,rawCorr_ch,'Color',[clr 0.75]*2/3,'LineWidth',1);
    ifsig = double(rawCorr_ch>=shuf_ch_high');
    ifsig(ifsig==0) = NaN;
    plot(T,rawCorr_ch.*ifsig,'Color',clr,'LineWidth',2)
    
    
    legend([a,b.mainLine],{'real',['shuffled x',num2str(niters)]},'box','off')
    xlabel('time from pics (ms)')
    ylabel(['correlation OFC:CN (',which_level,' PPD)'])
    
%     xlim([-500 1000])
    xlim([-1000 1000])
    
end


end


