 function bhvdata = load_cleanbhv(ses_info, dir_data, dir_bhv)
%load all bhv data from 'ses_info' table


% restrict to these variables:
varnames = {'experiment','session','subject','TrialNumber',...
    'jpg','lever','rt','ifreward','amnt','prob','trialtype',...
    'MLsacc','saccloc','subjval_expval','valbin_expval','location_target'};

% init bhvdata struct
bhvdata = struct();
for v = 1:length(varnames)
    bhvdata.(varnames{v}) = [];
end

nses = size(ses_info,1);
for s = 1:nses
    
    % get session details
    experiment = ses_info{s,1}{1};
    subject = ses_info{s,2}{1};
    num = ses_info{s,3}{1};
    date = strrep(ses_info{s,4}{1},'/','');
    
    session = strjoin({num,date},'_');
    
    % load session clean bhv file
    fileloc = fullfile(dir_data,experiment,subject,dir_bhv,session);
    %     fileloc = dir_bhv;
    date_fancy = [date(5:end),'-',date(1:2),'-',date(3:4)];
    
    if strcmp(subject,'George')
        temp = dir([fileloc,'/amnt*',date_fancy,'*_clean.mat']);
    elseif strcmp(subject,'Chap')
        temp = dir([fileloc,'/Experiment*',date_fancy,'*_clean.mat']);
    else
        error('wrong subject')
    end
    
    if length(temp)~=1
        error('check bhv files!!!')
    end
    filename = temp(1).name;
    
    load(fullfile(fileloc,filename));
    
    % update bhvdata
    ntr = length(trialinfo.TrialNumber);
    trialinfo.TrialNumber = reshape(trialinfo.TrialNumber,[],1); % transposed in some files
    
    trialinfo.MLsacc = trialinfo.MLsacc(:,1:2); % limit to 1st & 2nd saccade
    trialinfo.saccloc = trialinfo.saccloc(:,1:2); % limit to 1st & 2nd saccade
    
    trialinfo.experiment = repmat({experiment},ntr,1);
    trialinfo.subject = repmat({subject},ntr,1);
    trialinfo.session = repmat({session},ntr,1);
    
    
    for v = 1:length(varnames)
        bhvdata.(varnames{v}) = cat(1, bhvdata.(varnames{v}), trialinfo.(varnames{v}));
    end
    
end

% label correct trials = ch > unch by valbin
temp = bhvdata.valbin_expval;
temp(bhvdata.lever==1,:) = temp(bhvdata.lever==1,[2,1]);

bhvdata.chose_high = ones(size(temp,1),1);
bhvdata.chose_high(temp(:,1)==temp(:,2)) = 0;
bhvdata.chose_high(temp(:,1)<temp(:,2)) = -1;
bhvdata.chose_high(bhvdata.lever==0) = 0;


end

