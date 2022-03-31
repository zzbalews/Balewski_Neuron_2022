function [notch_filts,band_filts,bandnames] = make_filters(adfreq,path_figs)
% make notch & bandpass filters, assuming 'adfreq' is sampling freq (usually 1000Hz)

%% make notch filters: 3 @ 60Hz

Nfilts = 3; % number filters: 60, 120, 180
notch_filts = cell(1,Nfilts);

for f = 1:Nfilts
    
    freq = f*60;
    
    d = designfilt('bandstopfir', ... response type
        'FilterOrder',1000, ... filter order
        'CutoffFrequency1',freq-2, ... frequency constraints
        'CutoffFrequency2',freq+2, ...
        'SampleRate',adfreq); ... sampling rate
        
    notch_filts{f} = d;
    
end

% % show filters
% g=fvtool(notch_filts{:});
% xlim([0 250])
% print(strjoin({path_figs,['filt-00notch.png']},'/'),'-dpng')
% close(g)

%% make bandpass filters

freqbands = [... % all freq ranges
    2,4;... delta
    4,8;... theta
    8,12;... alpha
    12,30;... beta
    30,60;... gamma
    70,200];... high gamma
    
bandnames = {'01delta','02theta','03alpha','04beta','05gamma','06highgamma'};

Nfilts = size(freqbands,1);
band_filts = cell(1,Nfilts);

for f = 1:Nfilts
    
    d = designfilt('bandpassfir', ... response type
        'FilterOrder',1000, ... filter order
        'CutoffFrequency1', freqbands(f,1),... lower range
        'CutoffFrequency2', freqbands(f,2),... upper range
        'SampleRate', adfreq); ... sampling rate
        
    band_filts{f} = d;
    
end

% % show filters
% for f = 1:Nfilts
%     g=fvtool(band_filts{f});
%     xlim([0 250])
%     print(strjoin({path_figs,['filt-',bandnames{f},'.png']},'/'),'-dpng')
%     close(g)
% end

end

