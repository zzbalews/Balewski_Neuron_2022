function save_phase(outfile,pics_phase,choice_phase,sacc1_phase,sacc2_phase,...
    pics_range,choice_range,sacc1_range,sacc2_range)

save(outfile,'pics_phase','choice_phase','sacc1_phase','sacc2_phase',...
    'pics_range','choice_range','sacc1_range','sacc2_range','-v7.3');

end

