function show_LOO_acc(ses_info,LOO)

subjects = unique(ses_info.subject);
for s = 1:2
idx = strcmp(ses_info.subject,subjects{s});

figure; hold on

scatter(ones(sum(idx),1),LOO(idx),'filled',...
        'jitter','on',...
    'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',.4)

plot([0 2],[.25 .25],'k:')
ylim([0.15 .7])
xlim([.5 1.5])

ylabel({'Value decoder','leave-one-out accuracy'})
set(gca,'XTick',1,'XTickLabel',subjects{s},...
    'YTick',0:.25:1)

set(gcf,'Position',[100 100 150 300])
disp([subjects{s},' value decoding: LOO avg acc ',...
    num2str(mean(LOO(idx))),'+-',num2str(std(LOO(idx))/sqrt(sum(idx)))])

end


end

