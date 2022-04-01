function [dirdata,LOO] = load_dirdecoder(sync_point,ses_info, dir_data, dir_decode, bhvdata, get_LOO)


varnames = {'experiment','session','subject','TrialNumber','postprob'};

% init dirdata struct
dirdata = struct();
for v = 1:length(varnames)
    dirdata.(varnames{v}) = [];
end

LOO = struct();

nses = size(ses_info,1);
for s = 1:nses
    tic
    
    % get session details
    experiment = ses_info{s,1}{1};
    subject = ses_info{s,2}{1};
    num = ses_info{s,3}{1};
    date = strrep(ses_info{s,4}{1},'/','');
    
    session = strjoin({num,date},'_');
    
    % load session dir decoder output
    fileloc = fullfile(dir_data,experiment,subject,dir_decode,session);
    temp = dir([fileloc,'/*_',sync_point,'_freeLOO.mat']);
    filename = temp(1).name;
    
    load(fullfile(fileloc,filename),'alltrials_postprob_avg','t_mids');
    
    % update dirdata
    mini = struct();
    mini.postprob = permute(alltrials_postprob_avg,[1,2,4,3]);
    
    ntr = size(alltrials_postprob_avg,1);
    
    mini.experiment = repmat({experiment},ntr,1);
    mini.subject = repmat({subject},ntr,1);
    mini.session = repmat({session},ntr,1);
    mini.TrialNumber = (1:ntr)';
    
    for v = 1:length(varnames)
        dirdata.(varnames{v}) = cat(1, dirdata.(varnames{v}), mini.(varnames{v}));
    end
    
    % get session LOO accuracy
    if get_LOO
        load(fullfile(fileloc,filename),'LDA_output');
        
        tr = strcmp(LDA_output.trialgroup,'LOO');
        LOO.acc(s,:) = mean(mean(LDA_output.truelabel(tr,:,:)==LDA_output.guess(tr,:,:),3),1);
        LOO.subject(s,1) = {subject};
        
        % get LOO trial number & rt
        [~,t0] = min(abs(t_mids-0));
        LOO_trials = permute(LDA_output.trialnumber(tr,t0,:),[1,3,2]);
        
        rt = bhvdata.rt(strcmp(bhvdata.session,session));
        LOO_rt = rt(LOO_trials);
        
        % within each LOO iter, median split trials by 25% quantiles
        rt_splits = quantile(LOO_rt,[0,.25,.5,.75,1]);
        
        % get LOO accuracy by fast and slow
        [~,nmids,niters] = size(LDA_output.trialnumber);
        acc_quants = nan(niters,nmids,4);
       
       
        LDA_acc = LDA_output.truelabel(tr,:,:)==LDA_output.guess(tr,:,:);
        
        for iter = 1:niters
            
            for q = 1:4
                
                idx = LOO_rt(:,iter) > rt_splits(q,iter) & ...
                    LOO_rt(:,iter) <= rt_splits(q+1,iter);
                
                acc_quants(iter,:,q) = mean(LDA_acc(idx,:,iter));
            end
        end
       
        % update
        LOO.acc_quants(s,:,:) = mean(acc_quants);
        
    end
    toc
end

dirdata.t_mids = t_mids;
LOO.t_mids = t_mids;

% relabel by chosen direction
dirdata.postprob_ch = dirdata.postprob(:,:,1);
idxR = bhvdata.lever==1;
dirdata.postprob_ch(idxR,:) = dirdata.postprob(idxR,:,2);

end

