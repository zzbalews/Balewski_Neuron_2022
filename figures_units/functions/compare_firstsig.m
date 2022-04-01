function compare_firstsig(glmdata)


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
valcpd = glmdata.cpd(:,t_idx,valterm);
valsig = glmdata.sigunit(:,t_idx,valterm);
valunits = sum(valsig,2)>0;

% get first sig bin of value units
firstsig = nan(nunits,1);
peakcpd = nan(nunits,1);
for u = find(valunits)'
    
    temp = find(valsig(u,:),1);
    firstsig(u) = t_mids(temp + t_idx(1) - 1);
    
    peakcpd(u) = max(valcpd(u,:));
end

% relabel region
regionID = ones(nunits,1);
regionID(strcmp(glmdata.region,'CN')) = -1;

for s = 1:nsubj
    
    % restrict units
        subject = subject_names{s}
        
        idx =  strcmp(glmdata.subject,subject) & ...
            valunits;
        
        % terms
        Y = firstsig(idx);
        Y_norm = (Y - mean(Y))./std(Y);
        
        X = [regionID(idx),logit(peakcpd(idx))];
        X_norm = [ones(sum(idx),1), (X - mean(X))./std(X)];
        
        [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, Y_norm);
        B_hat
        coeff_p
        CPD
        
        
        
    
end






end

