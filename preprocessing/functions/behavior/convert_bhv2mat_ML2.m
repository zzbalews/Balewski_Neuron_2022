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
%   str_model_name = models to use for value fit, e.g. 'linear discount + first saccade'

%% skip if output exists!
%
% if exist(strjoin({path_out,outfile},'/'))
%     disp('processed bhv exists! skipping...')
%     return
% end

%% load raw bhv & extract trial info

lever_opts = [-1,1];

% init trial info struct
trialinfo = struct();

for i = 1:length(rawfiles)
    
    % init mini trial info struct
    smallinfo = struct();
    
    % filename
    BHVinput = strjoin({path_raw,rawfiles{i}},'/');
    
    % load raw bhv
    bhvdata = mlread(BHVinput);
    
    % trial #
    ntr = length(bhvdata);
    smallinfo.TrialNumber = (1:ntr)';
    
    % ML absoluate trial start
    smallinfo.MLAbsoluteTrialStartTime = [bhvdata.AbsoluteTrialStartTime]';
    
    % trial error
    smallinfo.MLTrialError = [bhvdata.TrialError]'; %0=correct; 1=no grab; 5=early; 2=late
    
    % ML trial start
    smallinfo.MLstart = pull_eventcode({bhvdata.BehavioralCodes},9);
    
    % ML trial stop
    smallinfo.MLstop = pull_eventcode({bhvdata.BehavioralCodes},18);
    
    % ML fix on
    smallinfo.MLfixon = pull_eventcode({bhvdata.BehavioralCodes},40);
    
    % ML pics on
    smallinfo.MLpicson = pull_eventcode({bhvdata.BehavioralCodes},20);
    
    % ML choice
    leftlever = pull_eventcode({bhvdata.BehavioralCodes},23);
    rightlever = pull_eventcode({bhvdata.BehavioralCodes},24);
    temp = nansum([leftlever,rightlever],2);
    temp(temp==0) = NaN;
    smallinfo.MLchoice = temp;
    
    % lever & rt
    smallinfo.lever = ~isnan(leftlever)*-1 + ~isnan(rightlever); % -1 = left, +1 = right
    smallinfo.rt = smallinfo.MLchoice - smallinfo.MLpicson; % time from pics on to choice
    
    % ML reward start
    temp = cellfun(@(x) x.StartTimes,{bhvdata.RewardRecord},'UniformOutput',false);
    temp(cellfun(@(x) isempty(x), temp)) = {NaN};
    smallinfo.MLrewardstart = [temp{:}]';
    
    % ML reward stop
    temp = cellfun(@(x) x.EndTimes,{bhvdata.RewardRecord},'UniformOutput',false);
    temp(cellfun(@(x) isempty(x), temp)) = {NaN};
    smallinfo.MLrewardstop = [temp{:}]';
    
    % if reward given
    smallinfo.ifreward = ~isnan(smallinfo.MLrewardstart); % 0/1 reward given
    
    % get left/right stim info
    temp = cellfun(@(x) x.CurrentConditionInfo(1).amnt,{bhvdata.TaskObject},'UniformOutput',false)';
    smallinfo.amnt = cell2mat(temp);
    
    temp = cellfun(@(x) x.CurrentConditionInfo(1).prob,{bhvdata.TaskObject},'UniformOutput',false)';
    smallinfo.prob = cell2mat(temp);
    
    % assign unique id for each picture
    smallinfo.jpg = get_jpg(smallinfo.amnt,smallinfo.prob); % jpg picture ID
    
    % trial type
    temp = cellfun(@(x) x.CurrentConditionInfo(1).trialtype,{bhvdata.TaskObject},'UniformOutput',false)';
    smallinfo.trialtype = strcmp(temp,'free')+1;
    
    
    %% from eye info, add first sacc after pics
    T = bhvdata(1).AnalogData.SampleInterval;
    
    eyes = cellfun(@(x) x.Eye, {bhvdata.AnalogData}, 'UniformOutput',false);
    eyes_diff = cellfun(@(x) [0 0;diff(x)], eyes, 'UniformOutput',false);
    eyes_velocity = cellfun(@(x) fastsmooth(sqrt(x(:,1).^2 + x(:,2).^2),3), eyes_diff,'UniformOutput',false);
    
    %     thresh_sacc = 0.4; % use this for Chap
    thresh_sacc = 0.65; % use this for George
    
    smallinfo.MLsacc = nan(ntr,2);
    smallinfo.saccloc = nan(ntr,2);
    smallinfo.saccN = nan(ntr,1);
    
    smallinfo.location_raw = nan(ntr,4001);
    smallinfo.location_target = nan(ntr,4001);
    
    
    % show all x/y position for all trials
    figure; hold on
    daspect([1 1 1])
    xlim([-20 20])
    ylim([-10 10])
    
    
    
    for tr = smallinfo.TrialNumber'
        % get saccade candidates
        sacc_start = find([-1; diff(eyes_velocity{tr} > thresh_sacc)]==1);
        sacc_stop = find([1; diff(eyes_velocity{tr} > thresh_sacc)]==-1);
        sacc_max = nan(size(sacc_start));
        
        for s = 1:length(sacc_start)
            [~,idx] = max(eyes_velocity{tr}(sacc_start(s):sacc_stop(s)));
            sacc_max(s) = idx + sacc_start(s) -1;
        end
        
        
        % find 1st 2 sacc after pics on
        idx = find((T*sacc_max) >= (smallinfo.MLpicson(tr)) & ...
            (T*sacc_max) < smallinfo.MLchoice(tr));
        
        if ~isempty(idx)
            
            % time of sacc
            sacc_time = reshape([T*sacc_max(idx) - smallinfo.MLpicson(tr)],[],1);
            
            % distance between points
            sacc_move = sqrt(sum((eyes{tr}(sacc_max(idx)+10,:) - eyes{tr}(sacc_max(idx)-10,:)).^2,2));
            
            % get sacc loc
            temp = eyes{tr}(:,1)';
            x_pos = (mean(temp((5:15)+sacc_max(idx)),2)>0)*2-1;
            
            % restrict to sacc > 80 apart & distance > 3
            keep = find([100; diff(sacc_time)]>80 & sacc_move > 3);
            smallinfo.saccN(tr) = length(keep);
            
            % save
            if length(keep)==1
                smallinfo.MLsacc(tr,1) = sacc_time(keep);
                smallinfo.saccloc(tr,1) = x_pos(keep);
            elseif length(keep)>1
                smallinfo.MLsacc(tr,1:2) = sacc_time(keep(1:2));
                smallinfo.saccloc(tr,1:2) = x_pos(keep(1:2));
            end
        end
        
        
        %
        %         if tr<10
        % %
        %             figure; hold on
        %             x = T*(1:length(eyes_velocity{tr}));
        %             plot(x,eyes_velocity{tr})
        %             plot(smallinfo.MLpicson(tr)*[1 1],[0 2])
        %             plot(x([1,end]),.25*[1 1])
        %             plot(x([1,end]),.5*[1 1])
        %             plot(x([1,end]),.75*[1 1])
        %             plot(x([1,end]),[1 1])
        %             ylim([0 5])
        %
        %             scatter(smallinfo.MLpicson(tr) + smallinfo.MLsacc_bybin{tr},zeros(size(smallinfo.MLsacc_bybin{tr})))
        %         end
        
        
        % save relevant eye position
        if ~isnan(smallinfo.MLpicson(tr))
            
            fix_x = mean(eyes{tr}(round(smallinfo.MLpicson(tr)/2) + (-300:-50),1));
            bounds_L = fix_x - [15 7];
            bounds_R = fix_x + [7 15];
            
            t = 2*(1:size(eyes{tr},1));
            t_shift = t-smallinfo.MLpicson(tr);
            
            raw_x_pos = interp1(t_shift,eyes{tr}(:,1),[-2000:2000]);
            raw_x_pos(abs(raw_x_pos)>15) = NaN;
            
            raw_y_pos = interp1(t_shift,eyes{tr}(:,2),[-2000:2000]);
%             raw_x_pos(abs(raw_y_pos)>4) = NaN;
            raw_x_pos(abs(raw_y_pos+4)>4) = NaN;
            
            smallinfo.location_raw(tr,:) = raw_x_pos;
            
%             in_left = raw_x_pos>=bounds_L(1) & raw_x_pos<=bounds_L(2) & ...
%                 raw_y_pos>=-3 & raw_y_pos<=3;
%             in_right = raw_x_pos>=bounds_R(1) & raw_x_pos<=bounds_R(2) & ...
%                 raw_y_pos>=-3 & raw_y_pos<=3;
            
            in_left = raw_x_pos>=bounds_L(1) & raw_x_pos<=bounds_L(2) & ...
                raw_y_pos+4>=-3 & raw_y_pos+4<=3;
            in_right = raw_x_pos>=bounds_R(1) & raw_x_pos<=bounds_R(2) & ...
                raw_y_pos+4>=-3 & raw_y_pos+4<=3;
            
            smallinfo.location_target(tr,:) = -in_left + in_right;
            
            temp = -in_left + in_right;
            idx = temp==-1;
            plot(raw_x_pos(idx),raw_y_pos(idx),'Color',[1 0 0 .01]);
            idx = temp==0;
            plot(raw_x_pos(idx),raw_y_pos(idx),'Color',[0 0 0 .01]);
            idx = temp==1;
            plot(raw_x_pos(idx),raw_y_pos(idx),'Color',[0 0 1 .01]);
            
        end
    end
    
    
    
    
    % combine multiple files
    varnames = fieldnames(smallinfo);
    for v = 1:length(varnames)
        if i==1
            trialinfo.(varnames{v}) = smallinfo.(varnames{v});
        else
            trialinfo.(varnames{v}) = cat(1,trialinfo.(varnames{v}),smallinfo.(varnames{v}));
        end
    end
    
end

%
% %% sanity check 1st & 2nd saccade
% figure;
%
% free_tr = find(trialinfo.trialtype==2 & trialinfo.lever~=0);
%
% subplot(1,4,1);
% trialinfo.saccN(isnan(trialinfo.saccN)) = 0;
% histogram(trialinfo.saccN(free_tr),'BinWidth',1,'FaceColor',[.5 .5 .5])
%
% xlim([0 10])
% xlabel('# of saccades/trial')
% ylabel('trial count')
%
% free_tr = find(trialinfo.trialtype==2 & trialinfo.lever~=0);
%
% subplot(1,4,2:3); hold on
% clrs = [0 .7 0; .7 0 .7];
%
% for i = 1:2
%
%     histogram(trialinfo.MLsacc(free_tr,i),'BinWidth',10,'FaceColor',clrs(i,:));
%
% end
%
% xlim([0 900])
% xlabel('saccade from pics on (ms)')
% ylabel('trial count')
%
% subplot(1,4,4); hold on
% for i = 1:2
%     temp = trialinfo.saccloc(free_tr,i) == trialinfo.lever(free_tr);
%     temp = temp(~isnan(trialinfo.saccloc(free_tr,i)));
%
%     bar(i,mean(temp),'FaceColor',clrs(i,:));
% end
% set(gca,'XTick',1:2,'XTickLabel',{'1st','2nd'})
% xlabel('saccade')
% ylabel({'proportion of trials','to lever direction'})
%
% set(gcf,'Position',[100 100 1200 500])

%%

% fix trial numbers if needed to combine multiple files
if length(rawfiles)>1
    ntr = length(trialinfo.TrialNumber);
    trialinfo.TrialNumber = 1:ntr;
end

%% if LFP, add pl2 timestamps
if ~isempty(lfpfile)
    
    % init within trial PL2 timing
    timepoints = {'start','stop','fixon','picson','choice'};%,'rewardstart','rewardstop'};
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
    
    % quick fix: in a couple of files, sequence [18;19] should be [55,19]
    temp = find(sv(1:end-1)==18);
    badidx = temp(find(sv(temp+1)==19));
    sv(badidx) = 55;
    
    % get trial start and stops
    starts = find(sv==9);
    stops = find(sv==18);
    
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
        temp = ts_tr(sv_tr==40);
        trialinfo.PL2fixon(tr) = temp(end);
        
        % PL2 pics on
        temp = ts_tr(sv_tr==20);
        if ~isempty(temp)
            trialinfo.PL2picson(tr) = temp;
        end
        
        % PL2 choice/lever
        temp = {ts_tr(sv_tr==23),ts_tr(sv_tr==24)};
        if sum(cellfun(@(x) ~isempty(x),temp))>0 % lever selected!
            
            trialinfo.PL2choice(tr) = [temp{:}];
            
            %             % PL2 reward
            %             if length(bhvdata(tr).RewardRecord.StartTimes)>0
            %                 trialinfo.PL2rewardstart(tr) = bhvdata(tr).RewardRecord.StartTimes(end);
            %                 trialinfo.PL2rewardstop(tr) = bhvdata(tr).RewardRecord.EndTimes(end);
            %             end
            
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

