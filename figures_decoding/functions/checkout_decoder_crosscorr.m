function checkout_decoder_crosscorr(bhvdata,decdata1,decdata2,keep_tr,which_level)
%%
% useful vars

subject_names = unique(bhvdata.subject);
nsubj = length(subject_names);
BB = {[-500,0],[0,500],[500,1000]};
figure;
for bnd = 1:3
    bounds = BB{bnd};%[-850,1200];
    
    t_mids = decdata1.t_mids;
    t_subset = find(t_mids>=bounds(1) & t_mids<bounds(2));
    T = t_mids(t_subset);
    T_crosscorr = -(length(T)-1):(length(T)-1);
    niters = 1000;
    
    pairs = [1,2; 1,3; 1,4;2,3;2,4;3,4];
    
    if strcmp(which_level,'chosen')
        clr = [1 0 0];
        ppd1 = 'ch_ppd';
        ppd2 = 'ch_ppd';
    elseif strcmp(which_level,'unchosen')
        clr = [0 0 1];
        ppd1 = 'unch_ppd';
        ppd2 = 'unch_ppd';
        
    elseif strcmp(which_level,'val+dir')
        clr = [0 .8 0];
        ppd1 = 'ch_ppd';
        ppd2 = 'postprob_ch';
    end
    
    
    for s = 1:nsubj
        subject = subject_names{s};
        idx = find(keep_tr & ...
            strcmp(bhvdata.subject, subject));
        
        ntr = length(idx);
        
        if ntr==0
            continue
        end
        
        OFC_PPD_ch = decdata1.(ppd1)(idx,t_subset);
        CN_PPD_ch = decdata2.(ppd2)(idx,t_subset);
        
        
        subplot(2,3,bnd+3*(s-1)); hold on
        
        title(subject)
        
        % show shuffled crosscorrelation distribution
        track_xcorr_shuf_avg = nan(niters,length(T_crosscorr));
        
        
        
        tic
        for i = 1:niters
            
            if mod(i,100)==0
                i
                toc
            end
            row_shuf = nan(ntr,1);
            for p = 1:size(pairs,1)
                idx_pair = find(min(bhvdata.valbin_expval(idx,:),[],2)==pairs(p,1) & ...
                    max(bhvdata.valbin_expval(idx,:),[],2)==pairs(p,2));
                row_shuf(idx_pair) = idx_pair(randperm(length(idx_pair)));
            end
            
            %         row_shuf = randperm(ntr);
            OFC_PPD_ch_shuf = OFC_PPD_ch(row_shuf,:);
            
            cnorm = normxcorr2(CN_PPD_ch,OFC_PPD_ch_shuf);
            track_xcorr_shuf_avg(i,:) = cnorm(ntr,:);
            
        end
        
        shuf_ch_low = quantile(track_xcorr_shuf_avg,.01,1);
        shuf_ch_med = quantile(track_xcorr_shuf_avg,.5,1);
        shuf_ch_high = quantile(track_xcorr_shuf_avg,.99,1);
        
        b=shadedErrorBar(T_crosscorr,shuf_ch_med,...
            [shuf_ch_med-shuf_ch_low; shuf_ch_high-shuf_ch_med],...
            'lineprops',{'Color',[.6 .6 .6]});
        
        
        % show real crosscorrelation
        
        cnorm = normxcorr2(CN_PPD_ch,OFC_PPD_ch);
        track_xcorr_avg = cnorm(ntr,:);
        
        
        plot(T_crosscorr,track_xcorr_avg,'Color',[1 0 0 .5],'LineWidth',1)
        ifsig = double(track_xcorr_avg > shuf_ch_high);
        ifsig(ifsig==0) = NaN;
        plot(T_crosscorr,track_xcorr_avg.*ifsig,'Color',[1 0 0],'LineWidth',2)
        YLIM = [-0.07 0.13];%YLIM=get(gca,'YLim');
        ylim(YLIM)
        plot([0 0],YLIM,'k:')
        set(gca,'XTick',-100:50:100,'XTickLabel',5*(-100:50:100),...
            'YTick',-0.05:0.05:.15);%linspace(YLIM(1),YLIM(2),5))
        xlabel('\leftarrow OFC leads                    lag (ms)                    CN leads \rightarrow')
        ylabel('cross correlation (normalized)')
        xlim([-50 50])
        %     xlim([-100 100])
    end
    
    
end
end


