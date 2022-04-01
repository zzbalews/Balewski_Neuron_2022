function show_dirdec_heatmap(sync_point,bhvdata,dirdata,keep_tr)

subject_names = unique(dirdata.subject);
nsubj = length(subject_names);

t_mids = dirdata.t_mids;
nmids = length(t_mids);

[~,t0] = min(abs(t_mids-0));
[~,tfix] = min(abs(t_mids+750));

if strcmp(sync_point,'pics')
[~,t_start] = min(abs(t_mids+500));
[~,t_stop] = min(abs(t_mids-1500));

elseif strcmp(sync_point,'choice')
  [~,t_start] = min(abs(t_mids+1500));
[~,t_stop] = min(abs(t_mids-500));  
    
end
% define colormap
X = redbluecmap;
X = min(X([11:-1:6,6:-1:1],:)+.2,1);
X_smooth = mat2colormap(X);

figure;
for S = 1:nsubj
    % restrict trials
    subject = subject_names{S};
    idx = find(strcmp(bhvdata.subject,subject_names{S}) & ...
        keep_tr);
    
    ntr = length(idx);
   
    % order by RT
    rt = bhvdata.rt(idx);
    [rt_sort,row_sort] = sort(rt);
    
    % plot
    subplot(1,2,S)
    imagesc(dirdata.postprob_ch(idx(row_sort),:))
%     imagesc(dirdata.postprob_ch_norm(idx(row_sort),:))
    hold on
    
    % add t0 line
    plot([t0 t0],0.5 + [0,ntr],'Color',[.5 .5 .5],'LineWidth',2)
    
    if strcmp(sync_point,'pics')
  
    plot([tfix tfix],0.5 + [0,ntr],'k:','LineWidth',2)   % add fix line
    xlabel('time from pics (s)')
    elseif strcmp(sync_point,'choice')
        rt_sort = -rt_sort;
      xlabel('time from choice (s)')  
    end
    
    % add RT line
    plot(t0 + rt_sort/5,1:ntr,'k','LineWidth',2)
    
    % format
    xlim([t_start,t_stop])
    caxis([0 1])
%     caxis([-4 4])
    temp = t0:-100:1;
    ticks = [temp(end:-1:2),t0:100:nmids];
    set(gca,'YTick',0:2500:ntr,...
        'XTick',ticks,'XTickLabel',(t_mids(ticks)+1)/1000,...
        'box','off')
    
    
    ylabel('trials')
%     title(subject)

% cc=colorbar;
% title(cc,{'chosen','post prob'})
% title(cc,{'chosen','post prob','\sigma fix'})
end

colormap(X_smooth)
% colormap([1 0 0; 1 1 1; 0 0 1])


end

