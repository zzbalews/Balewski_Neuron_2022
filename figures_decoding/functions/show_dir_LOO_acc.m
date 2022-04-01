function show_dir_LOO_acc(sync_point,LOO, bhvdata)

figure;

subjects = unique(LOO.subject);
for s = 1:2
    
   idx = strcmp(LOO.subject,subjects{s});
   
   subplot(1,2,s); hold on
   
   plot([-2000 2000],[.5 .5],'k:')
   plot([0 0],[0 1],'k:')
   
   plot(LOO.t_mids, LOO.acc(idx,:), 'Color', [.5 .5 .5 .3],'LineWidth',1)
   plot(LOO.t_mids, mean(LOO.acc(idx,:)), 'Color', [0 0 0], 'LineWidth',2)
   
   set(gca,'XTick',-2000:500:2000,'XTickLabel',-2:.5:2)
   xlim([-500, 1500])
   ylim([.4 1])
   xlabel('Time from pics (s)')
   ylabel('Direction deocder LOO accuracy')
   
    
end
end

