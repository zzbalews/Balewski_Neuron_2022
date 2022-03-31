function  working_signal = load_notched_LFP(fname,chname,notch_filts,path_figs,session_name)
% load LFP & apply notch filters

% load raw LFP
[adfreq,~,~,~,working_signal] = plx_ad(fname,chname);
% save_powerspec(working_signal,path_figs,[session_name,'_',chname],adfreq,'raw');

for f = 1:length(notch_filts)
    working_signal = fftfilt(notch_filts{f}.Coefficients,working_signal);
end
% save_powerspec(working_signal,path_figs,[session_name,'_',chname],adfreq,'notchfilt');


end

