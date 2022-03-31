function outfile = convert_bhv2mat_ML1(...
    path_raw,rawfiles,path_pl2,lfpfile,...
    path_out,outfile,str_model_names)
%convert_bhv2mat_ML1(path_raw,rawfiles,path_pl2,lfpfile,path_out,outfile)
%   pull trial info and eventcode timing
%   path_raw = location of *.bhv files
%   raw_files = names of *.bhv files
%   path_pl2 = location of LFP file
%   lfpfile = name of LFP file
%   path_out = location of processed trialinfo files
%   outfile = name of processed trialinfo file
%   session_name = recording session name
%   str_model_name = model to use for value fit, e.g. 'linear discount + first saccade'

%% skip if output exists!

% if exist(strjoin({path_out,outfile},'/'))
%     disp('processed bhv exists! skipping...')
%     return
% end

%% get stim info
load('task_params_Chap.mat')
lever_opts = [-1,1];

%% load raw bhv
% init bhv file struct
bhvdata = struct();
bhvdata_vars = {'AbsoluteTrialStartTime',...
    'ConditionNumber','TrialError','CodeNumbers','CodeTimes',...
    'RewardRecord','AnalogData'};
for v = 1:length(bhvdata_vars)
    bhvdata.(bhvdata_vars{v}) = [];
end

for i = 1:length(rawfiles) % combine bhv files, if multiple for session
    % filenames
    BHVinput = strjoin({path_raw,rawfiles{i}},'/');
    
    % load raw bhv
    temp = bhv_read(BHVinput);
    
    % combine files
    temp.CodeNumbers = temp.CodeNumbers';
    temp.CodeTimes = temp.CodeTimes';
    temp.RewardRecord = temp.RewardRecord';
    temp.AnalogData = temp.AnalogData';
    for v = 1:length(bhvdata_vars)
        bhvdata.(bhvdata_vars{v}) = cat(1,bhvdata.(bhvdata_vars{v}),temp.(bhvdata_vars{v}));
    end
    
    % keep TaskObject
    if i==1
        bhvdata.TaskObject = temp.TaskObject;
        bhvdata.AnalogInputFrequency = temp.AnalogInputFrequency;
    end
end

%% extract trial info
% init trial info struct
trialinfo = struct();

% trial #
ntr = length(bhvdata.AbsoluteTrialStartTime);
trialinfo.TrialNumber = (1:ntr)';

% ML trial start
trialinfo.MLAbsoluteTrialStartTime = bhvdata.AbsoluteTrialStartTime;

% trial error
trialinfo.MLTrialError = bhvdata.TrialError;

% init within trial ML timing
timepoints = {'start','stop','fixon','picson','choice','rewardstart','rewardstop'};
for i = 1:length(timepoints)
    trialinfo.(['ML',timepoints{i}]) = nan(ntr,1);
end
% trialinfo.MLstart = nan(ntr,1);
% trialinfo.MLstop = nan(ntr,1);
% trialinfo.MLfixon = nan(ntr,1);
% trialinfo.MLpicson = nan(ntr,1);
% trialinfo.MLchoice = nan(ntr,1);
% trialinfo.MLrewardstart = nan(ntr,1);
% trialinfo.MLrewardstop = nan(ntr,1);

% init stim info (left, right)
trialinfo.jpg = nan(ntr,2); % jpg picture ID

% init choice info
trialinfo.lever = zeros(ntr,1); % -1 = left, +1 = right
trialinfo.rt = nan(ntr,1); % time from pics on to choice

% init outcome info
trialinfo.ifreward = nan(ntr,1); % 0/1 reward given

% pull event codes and times within each trial and populate fields
for tr = trialinfo.TrialNumber'
    
    % pull events; remove repeated codes (sent x3 for backup)
    allcodes = bhvdata.CodeNumbers{tr};
    keep = [1; diff(allcodes)~=0];
    CodeNumbers = allcodes(keep==1); % event types
    CodeTimes = bhvdata.CodeTimes{tr}(keep==1); % corresponding timestamps
    
    % ML start trial
    trialinfo.MLstart(tr) = CodeTimes(CodeNumbers==9);
    
    % ML stop trial
    trialinfo.MLstop(tr) = CodeTimes(CodeNumbers==18);
    
    % ML fix on
    temp = CodeTimes(CodeNumbers==20);
    trialinfo.MLfixon(tr) = temp(end);
    
    % ML pics on
    temp = CodeTimes(CodeNumbers==4);
    if ~isempty(temp)
        trialinfo.MLpicson(tr) = temp;
    end
    
    % ML choice/lever
    temp = {CodeTimes(CodeNumbers==15),CodeTimes(CodeNumbers==10)};
    if sum(cellfun(@(x) ~isempty(x),temp))>0 % lever selected!
        trialinfo.MLchoice(tr) = [temp{:}];
        trialinfo.lever(tr) = lever_opts(find(cellfun(@(x) ~isempty(x),temp)));
        trialinfo.rt(tr) = trialinfo.MLchoice(tr) - trialinfo.MLpicson(tr);
        
        % ML reward
        trialinfo.ifreward(tr) = sum(CodeNumbers==30);
        trialinfo.MLrewardstart(tr) = bhvdata.RewardRecord(tr).RewardOnTime(end);
        trialinfo.MLrewardstop(tr) = bhvdata.RewardRecord(tr).RewardOffTime(end);
        
    end
    
    % jpg ID and stim location
    cond = bhvdata.ConditionNumber(tr);
    stim1 = strsplit(bhvdata.TaskObject{cond,3},',');
    stim2 = strsplit(bhvdata.TaskObject{cond,4},',');
    
    temp = cellfun(@(x) strrep(x,'Sqr(1','Pic(stim17'), {stim1{1},stim2{1}},'UniformOutput',false);
    jpgs = cellfun(@(x) str2num(strrep(x,'Pic(stim','')),temp);
    
    locs = cellfun(@(x) str2num(strrep(x,')','')),{stim1{end},stim2{end}});
    if locs(1) < locs(2) % stim1 = right
        trialinfo.jpg(tr,:) = jpgs([2,1]);
    else % stim1 = left
        trialinfo.jpg(tr,:) = jpgs;
    end
    
end

% stim amount & probability
for i = 1:2
    trialinfo.amnt(:,i) = stims(trialinfo.jpg(:,i),1); % juice amount (pump V)
    trialinfo.prob(:,i) = stims(trialinfo.jpg(:,i),2); % prob of reward
end

% trial type
trialinfo.trialtype = sum(trialinfo.jpg~=17,2);

%% from eye info, add first sacc after pics
freq = bhvdata.AnalogInputFrequency;
T = 1000/freq;

eyes = cellfun(@(x) x.EyeSignal, bhvdata.AnalogData, 'UniformOutput',false);
eyes_diff = cellfun(@(x) [0 0;diff(x)], eyes, 'UniformOutput',false);
eyes_velocity = cellfun(@(x) fastsmooth(sqrt(x(:,1).^2 + x(:,2).^2),3), eyes_diff,'UniformOutput',false);

thresh_sacc = 0.4; % use this for Chap
%     thresh_sacc = 0.65; % use this for George


trialinfo.MLsacc = nan(ntr,5);
trialinfo.saccloc = nan(ntr,5);
trialinfo.saccN = nan(ntr,1);

trialinfo.location_raw = nan(ntr,4001);
trialinfo.location_target = nan(ntr,4001);


% show all x/y position for all trials
figure; hold on
daspect([1 1 1])
xlim([-20 20])
ylim([-10 10])


for tr = trialinfo.TrialNumber'
    
    sacc_start = find([-1; diff(eyes_velocity{tr} > thresh_sacc)]==1);
    sacc_stop = find([1; diff(eyes_velocity{tr} > thresh_sacc)]==-1);
    sacc_max = nan(size(sacc_start));
    for s = 1:length(sacc_start)
        [~,idx] = max(eyes_velocity{tr}(sacc_start(s):sacc_stop(s)));
        sacc_max(s) = idx + sacc_start(s) -1;
    end
    
    
    % find 1st 2 sacc after pics on
    idx = find((T*sacc_max) >= (trialinfo.MLpicson(tr)) & ...
        (T*sacc_max) < trialinfo.MLchoice(tr));
    
    if ~isempty(idx)
        
        % time of sacc
        sacc_time = reshape([T*sacc_max(idx) - trialinfo.MLpicson(tr)],[],1);
        
        % distance between points
        sacc_move = sqrt(sum((eyes{tr}(sacc_max(idx)+10,:) - eyes{tr}(sacc_max(idx)-10,:)).^2,2));
        
        % get sacc loc
        temp = eyes{tr}(:,2)';
        y_pos = (mean(temp((5:15)+sacc_max(idx)),2)<0)*2-1;
        
        % restrict to sacc > 80 apart & distance > 3
        keep = find([100; diff(sacc_time)]>80 & sacc_move > 3);
        trialinfo.saccN(tr) = length(keep);
        
        % save
        if length(keep)==1
            trialinfo.MLsacc(tr,1) = sacc_time(keep);
            trialinfo.saccloc(tr,1) = y_pos(keep);
        elseif length(keep)>1
            trialinfo.MLsacc(tr,1:2) = sacc_time(keep(1:2));
            trialinfo.saccloc(tr,1:2) = y_pos(keep(1:2));
        end
    end
    
    %     if tr<10
    %
    %         figure; hold on
    %         x = T*(1:length(eyes_velocity{tr}));
    %         plot(x,eyes_velocity{tr})
    %         ylim([0 4])
    %         plot(trialinfo.MLpicson(tr)*[1 1],[0 2])
    %         plot(x([1,end]),.25*[1 1])
    %         plot(x([1,end]),.5*[1 1])
    %         plot(x([1,end]),.75*[1 1])
    %         plot(x([1,end]),[1 1])
    %     end
    
    % save relevant eye position
    if ~isnan(trialinfo.MLpicson(tr))
        
        fix_y = mean(eyes{tr}(round(trialinfo.MLpicson(tr)/2) + (-300:-50),2));
        bounds_L = fix_y + [2 8];
        bounds_R = fix_y - [8 2];
        
        t = 2*(1:size(eyes{tr},1));
        t_shift = t-trialinfo.MLpicson(tr);
        
        raw_y_pos = interp1(t_shift,eyes{tr}(:,2),[-2000:2000]);
        raw_y_pos(abs(raw_y_pos)>15) = NaN;
        
        raw_x_pos = interp1(t_shift,eyes{tr}(:,1),[-2000:2000]);
        raw_y_pos(abs(raw_x_pos)>4) = NaN;
        
        
        trialinfo.location_raw(tr,:) = raw_y_pos;
        
        in_left = raw_x_pos>=-3 & raw_x_pos<=3 & ...
            raw_y_pos>=bounds_L(1) & raw_y_pos<=bounds_L(2);
        in_right = raw_x_pos>=-3 & raw_x_pos<=3 & ...
            raw_y_pos>=bounds_R(1) & raw_y_pos<=bounds_R(2);
        
        trialinfo.location_target(tr,:) = -in_left + in_right;
        
        temp = -in_left + in_right;
        idx = temp==-1;
        plot(raw_x_pos(idx),raw_y_pos(idx),'Color',[1 0 0 .01]);
        idx = temp==0;
        plot(raw_x_pos(idx),raw_y_pos(idx),'Color',[0 0 0 .01]);
        idx = temp==1;
        plot(raw_x_pos(idx),raw_y_pos(idx),'Color',[0 0 1 .01]);
        
        
    end
    
    
end
% %%
% figure;
% for i = 1:5
% subplot(5,1,i);
%     histogram(trialinfo.MLsacc(:,i),'BinWidth',10)
%     title({['saccade #',num2str(i)],...
%         [num2str(mean(~isnan(trialinfo.MLsacc(trialinfo.lever~=0,i)))),' trials w/ sacc']})
%     xlim([0 1000])
%
% end


%% sanity check 1st & 2nd saccade
figure;

free_tr = find(trialinfo.trialtype==2 & trialinfo.lever~=0);

subplot(1,4,1);
trialinfo.saccN(isnan(trialinfo.saccN)) = 0;
histogram(trialinfo.saccN(free_tr),'BinWidth',1,'FaceColor',[.5 .5 .5])

xlim([0 10])
xlabel('# of saccades/trial')
ylabel('trial count')

free_tr = find(trialinfo.trialtype==2 & trialinfo.lever~=0);

subplot(1,4,2:3); hold on
clrs = [0 .7 0; .7 0 .7];

for i = 1:2
    
    histogram(trialinfo.MLsacc(free_tr,i),'BinWidth',10,'FaceColor',clrs(i,:));
    
end

xlim([0 900])
xlabel('saccade from pics on (ms)')
ylabel('trial count')

subplot(1,4,4); hold on
for i = 1:2
    temp = trialinfo.saccloc(free_tr,i) == trialinfo.lever(free_tr);
    temp = temp(~isnan(trialinfo.saccloc(free_tr,i)));
    
    bar(i,mean(temp),'FaceColor',clrs(i,:));
end
set(gca,'XTick',1:2,'XTickLabel',{'1st','2nd'})
xlabel('saccade')
ylabel({'proportion of trials','to lever direction'})

set(gcf,'Position',[100 100 1200 500])


%% if LFP, add pl2 timestamps
if ~isempty(lfpfile)
    
    % init within trial PL2 timing
    timepoints = {'start','stop','fixon','picson','choice','rewardstart','rewardstop','leverdir'};
    for i = 1:length(timepoints)
        trialinfo.(['PL2',timepoints{i}]) = nan(ntr,1);
    end
    
    % load LFP codes and timestamps
    [n,ts,sv] = plx_event_ts(strjoin({path_pl2,lfpfile},'/'), 'Strobed');
    
    % remove triplet codes, keeping first (like bhv)
    keep = [1; diff(sv)~=0];
    n = sum(keep); % # of events
    sv = sv(keep==1); % event codes
    ts = ts(keep==1); % timestamps, in sec
    
    % get trial start and stops
    starts = find(sv==9);
    stops = find(sv==18);
    
% % %         % hacky fix for Chap00_rec07_11172917: remove "trial 883" from raw pl2
% % %         % (this trial is missing between the two bhv files for this session)
% % %         starts = starts([1:881,883:end]);
% % %         stops = stops([1:881,883:end]);
    
    % get # trials
    ntrPL2 = length(stops);
    
    % same # trials?
    if ntr~=ntrPL2 % does not match bhv
        error('error: PL2 & ML trial #s do not match')
        return
    end
    
    % pull event codes and times within each trial and populate fields
    for tr = 1:ntrPL2
        
        % extract trial events and timestamps
        idx = starts(tr):stops(tr);
        sv_tr = sv(idx);
        ts_tr = ts(idx)*1000; % convert time to ms
        
        % PL2 start trial
        trialinfo.PL2start(tr) = ts_tr(sv_tr==9);
        
        % PL2 stop trial
        trialinfo.PL2stop(tr) = ts_tr(sv_tr==18);
        
        % PL2 fix on
        temp = ts_tr(sv_tr==20);
        trialinfo.PL2fixon(tr) = temp(end);
        
        % PL2 pics on
        temp = ts_tr(sv_tr==4);
        if ~isempty(temp)
            trialinfo.PL2picson(tr) = temp;
        end
        
        % PL2 choice/lever
        temp = {ts_tr(sv_tr==15),ts_tr(sv_tr==10)};
        if sum(cellfun(@(x) ~isempty(x),temp))>0 % lever selected!
            trialinfo.PL2choice(tr) = [temp{:}];
            trialinfo.PL2leverdir(tr) = find(cellfun(@(x) ~isempty(x), temp))*2-3;
            
            %             % PL2 reward
            %             trialinfo.PL2rewardstart(tr) = bhvdata.RewardRecord(tr).RewardOnTime(end);
            %             trialinfo.PL2rewardstop(tr) = bhvdata.RewardRecord(tr).RewardOffTime(end);
            %
        end
        
    end
    
    % omitted same trials?
    same_omit = isnan(trialinfo.MLchoice) == isnan(trialinfo.PL2choice);
    if sum(same_omit~=1)>0
        error('error: PL2 & ML skipped trials do not match')
    end
    
    % calculate saccade pl2times
    trialinfo.PL2sacc = trialinfo.PL2picson + trialinfo.MLsacc;
    
end

%% add psychometric fit

% use only completed free trials
keep = trialinfo.trialtype==2 & trialinfo.lever~=0 & ~isnan(trialinfo.saccloc(:,1));

% choice direction
lever = trialinfo.lever(keep); % -1=left, +1=right

% trial info (column 1 = left, 2 = right)
jpg = trialinfo.jpg(keep,:); % picture ID (17=empty space)
amnt = trialinfo.amnt(keep,:); % juice amount, in pump Voltage
prob = trialinfo.prob(keep,:); % prob of juice reward
firstsacc = trialinfo.saccloc(keep,1); % -1=left, +1=right


% do fit, according to model
forced_emptyside = trialinfo.jpg==17; % empty spaces on forced
for s = 1:length(str_model_names)
    if strcmp(str_model_names{s},'linear discount + first saccade')
        [w_names,w_fit,BIC,R2] = fit_discount_sacc('lin',lever,jpg,amnt,prob,firstsacc);
        subjval_lin = lindiscount(trialinfo.amnt,trialinfo.prob,w_fit(1));
        
        % add subjective value to trialinfo
        subjval_lin(forced_emptyside) = NaN;
        trialinfo.subjval_lin = subjval_lin;
        trialinfo.valbin_lin = bin_val(subjval_lin,4);
        
        
    elseif strcmp(str_model_names{s},'expected value + first saccade')
        [w_names,w_fit,BIC,R2] = fit_discount_sacc('expval',lever,jpg,amnt,prob,firstsacc);
        subjval_expval = expvaldiscount(trialinfo.amnt,trialinfo.prob);
        
        % add subjective value to trialinfo
        subjval_expval(forced_emptyside) = NaN;
        trialinfo.subjval_expval = subjval_expval;
        trialinfo.valbin_expval = bin_val(subjval_expval,4);
        
    else
        disp('ERROR: psychometric model name not recognized')
        return
    end
    
end

%% save file!

save(strjoin({path_out,outfile},'/'),'trialinfo');

end

