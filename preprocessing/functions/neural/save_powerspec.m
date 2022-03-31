function save_powerspec(working_signal,path_figs,chID,adfreq,stage)

f=figure;
periodogram(working_signal,[],[],adfreq);
xlim([0 250])
ylim([-100 60])

print(strjoin({path_figs,[chID,'_',stage,'.png']},'/'),'-dpng')
close(f);
end

