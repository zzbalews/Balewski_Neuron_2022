function save_mag(outfile,pics_mag,choice_mag,sacc1_mag,sacc2_mag,...
    pics_range,choice_range,sacc1_range,sacc2_range)

save(outfile,'pics_mag','choice_mag','sacc1_mag','sacc2_mag',...
    'pics_range','choice_range','sacc1_range','sacc2_range','-v7.3');

end

