function [mag, phs] = after_bandpass_LFP(working_signal, coeffs, lag)
% bandpass filter working signal, and account for timing lag; 
% smooth & zscore

% filter signal
band_signal = fftfilt(coeffs,working_signal);

% get analytic amplitude & phase (hilbert transform)
band_hilbert_mag = abs(hilbert(band_signal));
band_hilbert_phase = angle(hilbert(band_signal));

% adjust for time lag
mag = band_hilbert_mag(lag:end);
phs = band_hilbert_phase(lag:end);

% smooth magnitude by 50ms boxcar
mag = fastsmooth(mag,49,1,1);

% zscore magnitude
mu = mean(mag);
sig = std(mag);
mag = (mag-mu)/sig;

% add NaN at end of signal for incomplete trials
mag(end+1) = NaN;
phs(end+1) = NaN;

end

