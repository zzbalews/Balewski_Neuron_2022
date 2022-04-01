function [valdata,LOO] = load_valdecoder(ses_info, bhvdata, dir_data, dir_decode, get_LOO)

% # sessions
nses = size(ses_info,1);

% init valdata struct

varnames = {'experiment','session','subject','TrialNumber',...
    'ch_ppd','unch_ppd','na_ppd','ch_state','unch_state','na_state'};
valdata = struct();
for v = 1:length(varnames)
    valdata.(varnames{v}) = [];
end

% init LOO array
LOO = nan(nses,1);
for s = 1:nses
    tic
    
    % get session details
    experiment = ses_info{s,1}{1};
    subject = ses_info{s,2}{1};
    num = ses_info{s,3}{1};
    date = strrep(ses_info{s,4}{1},'/','');
    
    session = strjoin({num,date},'_');
    
    ntr = sum(strcmp(bhvdata.session,session));
    
    % load session dir decoder output
    fileloc = fullfile(dir_data,experiment,subject,dir_decode,session);
    temp = dir([fileloc,'/*_value_free_w400.mat']);
    filename = temp(1).name;
    
    load(fullfile(fileloc,filename),'ch_ppd','unch_ppd','na_ppd',...
        'ch_state','unch_state','na_state','t_mids','free_tr');
    
    nmids = length(t_mids);
    
    % update
    mini = struct();
    
    mini.ch_ppd = nan(ntr,nmids);
    mini.ch_ppd(free_tr,:) = ch_ppd;
    
    mini.unch_ppd = nan(ntr,nmids);
    mini.unch_ppd(free_tr,:) = unch_ppd;
    
    mini.na_ppd = nan(ntr,nmids,3);
    mini.na_ppd(free_tr,:,:) = na_ppd;
    
    mini.ch_state = nan(ntr,nmids);
    mini.ch_state(free_tr,:) = ch_state;
    
    mini.unch_state = nan(ntr,nmids);
    mini.unch_state(free_tr,:) = unch_state;
    
    mini.na_state = nan(ntr,nmids,3);
    mini.na_state(free_tr,:,:) = na_state;
    
    mini.experiment = repmat({experiment},ntr,1);
    mini.subject = repmat({subject},ntr,1);
    mini.session = repmat({session},ntr,1);
    mini.TrialNumber = (1:ntr)';
    
    for v = 1:length(varnames)
        valdata.(varnames{v}) = cat(1, valdata.(varnames{v}), mini.(varnames{v}));
    end
    
    valdata.t_mids = t_mids;
    % load LOO accuracy if requested
    if get_LOO
        
        temp = dir([fileloc,'/*_value_forcedLOO.mat']);
        filename = temp(1).name;
        
        load(fullfile(fileloc,filename),'m','t_mids');
        LOO(s) = max(m);
        
    end
    toc
end




end

