function compare_LOO_Nunits(glmdata, dir_data, ...
    dir_decode_valueOFC, dir_decode_valueCN, dir_decode_dirCN)

%%
region_names = unique(glmdata.region);
subject_names = unique(glmdata.subject);

nreg = length(region_names);
nsubj = length(subject_names);

t_mids = glmdata.t_mids;
nmids = length(t_mids);

t_idx = find(t_mids>=0 & t_mids<=500);

nunits = length(glmdata.unit_names);

% find all value units
valterm = 4;
valsig = glmdata.sigunit(:,t_idx,valterm);
valunits = sum(valsig,2)>0;

% find all direction nunits
dirterm = 3;
dirsig = glmdata.sigunit(:,t_idx,dirterm);
dirunits = sum(dirsig,2)>0;

valdec = struct();
dirdec = struct();

valdec.subject = cell(1,0);
dirdec.subject = cell(1,0);

valdec.nUnits = nan(1,0);
dirdec.nUnits = nan(1,0);

valdec.region = cell(1,0);
diredec.region = cell(1,0);

valdec.acc = nan(1,0);
dirdec.acc = nan(1,0);

j = 0;
k = 0;

for s = 1:nsubj
    
    subject = subject_names{s};
    idx_subj = strcmp(glmdata.subject,subject);
    session_names = unique(glmdata.session(idx_subj));
    nses = length(session_names);
    
    for ses = 1:nses
        % OFC value
        j = j+1;
        
        idx = idx_subj & valunits & ...
            strcmp(glmdata.region,'OFC') & ...
            strcmp(glmdata.session,session_names{ses});
        
        floc = fullfile(dir_data,subject,dir_decode_valueOFC,session_names{ses});
        temp = dir([floc,'/*LOOacc_w400.mat']);
        if ~isempty(temp)
            load(fullfile(floc,temp.name),'track_guess','track_label');
            m = reshape(mean(track_guess==track_label,3),1,[]);
            
            valdec.subject{j} = subject;
            valdec.region{j} = 'OFC';
            valdec.nUnits(j) = sum(idx);
            valdec.acc(j) = mean(m);
        end
        
        % CN value
        j = j+1;
        
        idx = idx_subj & valunits & ...
            strcmp(glmdata.region,'CN') & ...
            strcmp(glmdata.session,session_names{ses});
        
        floc = fullfile(dir_data,subject,dir_decode_valueCN,session_names{ses});
        temp = dir([floc,'/*LOOacc_w400.mat']);
        if ~isempty(temp)
            load(fullfile(floc,temp.name),'track_guess','track_label');
            m = reshape(mean(track_guess==track_label,3),1,[]);
            
            valdec.subject{j} = subject;
            valdec.region{j} = 'CN';
            valdec.nUnits(j) = sum(idx);
            valdec.acc(j) = mean(m);
        end
        
        % CN direction
        k = k+1;
        
        idx = idx_subj & dirunits & ...
            strcmp(glmdata.region,'CN') & ...
            strcmp(glmdata.session,session_names{ses});
        
        floc = fullfile(dir_data,subject,dir_decode_dirCN,session_names{ses});
        temp = dir([floc,'/*freeLOO.mat']);
        if ~isempty(temp)
            load(fullfile(floc,temp.name),'LDA_output');
            tr = strcmp(LDA_output.trialgroup,'LOO');
            m = mean(mean(LDA_output.truelabel(tr,:,:)==LDA_output.guess(tr,:,:),3),1);
            
            dirdec.subject{k} = subject;
            dirdec.region{k} = 'CN';
            dirdec.nUnits(k) = sum(idx);
            dirdec.acc(k) = max(m);
        end
    end
    
end

keep_rows = cellfun(@(x) ~isempty(x),valdec.subject);
valdec.subject = valdec.subject(keep_rows);
valdec.region = valdec.region(keep_rows);
valdec.nUnits = valdec.nUnits(keep_rows);
valdec.acc = valdec.acc(keep_rows);

%% plot
figure;

for s = 1:2
    
    subject = subject_names{s};
    
    % OFC value
    subplot(1,3,1); hold on
    idx = strcmp(valdec.subject,subject) & ...
        strcmp(valdec.region,'OFC');
    scatter(valdec.nUnits(idx), valdec.acc(idx),'filled')
    xlabel('# value OFC units')
    
    % CN value
    subplot(1,3,2); hold on
    idx = strcmp(valdec.subject,subject) & ...
        strcmp(valdec.region,'CN');
    scatter(valdec.nUnits(idx), valdec.acc(idx),'filled')
    xlabel('# value CN units')
    
    % CN direction
    subplot(1,3,3); hold on
    idx = strcmp(dirdec.subject,subject);
    scatter(dirdec.nUnits(idx), dirdec.acc(idx),'filled')
    xlabel('# direction CN units')
    
end

XLIM = [0 120];

for i = 1:3
    subplot(1,3,i);
    
    if i<3
        plot(XLIM,[.25 .25],'k:')
        ylim([.2 .75])
        ylabel('value decoder LOO accuracy')
    else
        plot(XLIM,[.5 .5],'k:')
        ylim([.45 1])
        ylabel('direction decoder LOO accuracy')
    end
    
    xlim(XLIM)
    
    if i==1
        legend(subject_names,'box','off','Location','northwest')
    end
        
end

Y = logit(valdec.acc');
N = length(Y);

temp = ones(N,1);
temp(strcmp(valdec.region,'CN')) = -1;
X = [temp, valdec.nUnits', temp.*valdec.nUnits'];

Y_norm = (Y-mean(Y))./std(Y);
X_norm = [ones(sum(keep_rows),1), (X-mean(X))./std(X)];

[ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, Y_norm);



end

