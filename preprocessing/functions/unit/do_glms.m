function outfname = do_glms(pl2dir, spkfname, minidata, units, outdir, which_val_model, which_trials)

%% determine output file, skip if exists
temp = strsplit(spkfname,'_');

% get hemisphere
hemisphere = temp{end}(1);

% generate out file
part1 = strjoin(temp([1:3]),'_');
part2 = strjoin({which_val_model, which_trials},'_');
part3 = strjoin(temp([8:9]),'_');
fname = strjoin({part1,part2,part3},'_');
outfname = strjoin({outdir, fname},'/');

if exist(outfname) % quit if files exist
    disp('glms done!')
    return
end
%% switch lever direction to ipsi/contra

if strcmp(hemisphere,'L')
    hemi_convert = 1;
elseif strcmp(hemisphere,'R')
    hemi_convert = -1;
end


%% load spk file
load(strjoin({pl2dir, spkfname},'/'),'SPKfr_units','SPKfr_names','t_mids','meta')
[~,nmids,nunits] = size(SPKfr_units);

%% format bhv as glm input

if strcmp(which_val_model,'expval')
    maxsubjvalue = max(minidata.subjval,[],2);
    diffsubjvalue = abs(diff(minidata.subjval,[],2));
elseif strcmp(which_val_model,'expval_bin')
    maxsubjvalue = max(minidata.valbin,[],2);
    diffsubjvalue = abs(diff(minidata.valbin,[],2));
else
    error('value model not recognized')
end


% use selected trial types
if strcmp(which_trials,'alltrials')
    varnames = {'constant','trialtype','leverdir','maxsubjvalue','trialnum'};
    
    predictors = table(...
        minidata.trialtype*2 - 3,... % forced = -1, free = +1
        minidata.lever * hemi_convert, ... % ipsi = -1, contra = +1
        maxsubjvalue,... % max subj value
        reshape(minidata.TrialNumber,[],1),... % trial number
        'VariableNames',varnames(2:end));
    
    % restrict to completed trials
    keep_tr = minidata.lever~=0;
    
elseif strcmp(which_trials,'free')
    varnames = {'constant','leverdir','maxsubjvalue','diffsubjvalue','trialnum'};
    
    predictors = table(...
        minidata.lever * hemi_convert, ... % ipsi = -1, contra = +1
        maxsubjvalue,... % max subj value
        diffsubjvalue,... % diff subj value
        reshape(minidata.TrialNumber,[],1),... % trial number
        'VariableNames',varnames(2:end));
    
    % restrict to completed free trials
    keep_tr = minidata.lever~=0 & minidata.trialtype==2;
    
elseif strcmp(which_trials,'freesimple')
    varnames = {'constant','leverdir','maxsubjvalue','trialnum'};
    
    predictors = table(...
        minidata.lever * hemi_convert, ... % ipsi = -1, contra = +1
        maxsubjvalue,... % max subj value
        reshape(minidata.TrialNumber,[],1),... % trial number
        'VariableNames',varnames(2:end));
    
   % restrict to completed free trials
    keep_tr = minidata.lever~=0 & minidata.trialtype==2;
    
    
else
    error('trial type invalid')
end

nterms = length(varnames);


% restrict trials
SPKfr_units = SPKfr_units(keep_tr,:,:);
predictors = predictors(keep_tr,:);

ntr = size(predictors,1);


% normalize
m = repmat(mean(SPKfr_units,1),ntr,1,1);
s = repmat(std(SPKfr_units,[],1),ntr,1,1);
SPKfr_units_norm = (SPKfr_units - m)./s;

predictors{:,:} = (predictors{:,:} - mean(predictors{:,:}))./std(predictors{:,:});

% make glm input
X = [ones(ntr,1),predictors{:,:}];

% get vif
vif = max(diag(inv(corrcoef(predictors{:,:}))));
vif
% median rt for session
medrt = median(minidata.rt(keep_tr));

%% set up for glms
track_beta = nan(nunits,nmids,nterms);
track_beta_p = nan(nunits,nmids,nterms);
track_model_p = nan(nunits,nmids);
track_cpd = nan(nunits,nmids,nterms);
track_expvar = nan(nunits,nmids);

%% do glms!
parfor T = 1:nmids
    for U = 1:nunits
        Y = SPKfr_units_norm(:,T,U);
        [B_hat, ~, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X,Y);
        track_beta(U,T,:) = B_hat;
        track_beta_p(U,T,:) = coeff_p;
        track_model_p(U,T) = model_p;
        track_cpd(U,T,:) = CPD;
        track_expvar(U,T) = (SSR_const - SSR)/SSR_const;
    end
end

%% clean up

output = struct();
temp = strsplit(spkfname,'_');
output.session = strjoin(temp(1:3),'_');
output.unit_names = SPKfr_names;
output.varnames = varnames;

temp = cellfun(@(x) find(strcmp(x,units.id)),output.unit_names);
output.unit_fr = units.fr(temp);

output.vif = vif;
output.beta = track_beta;
output.beta_p = track_beta_p;
output.model_p = track_model_p;
output.cpd = track_cpd;
output.expvar = track_expvar;

%% find sig units

thresh = 0.01; % min p value for beta & model
N = 4; % # consecutive time bins

sig_candidate = output.beta_p < thresh & output.model_p < thresh;
output.sigunit = nan(size(sig_candidate));

for U = 1:nunits
    for R = 1:nterms
        
        temp = sig_candidate(U,:,R);
        output.sigunit(U,:,R) = give_consec_seg(temp,N);
        
    end
end


%% save
save(outfname,'output','t_mids','meta','-v7.3')

end

