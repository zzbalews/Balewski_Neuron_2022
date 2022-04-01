function [example_units] = pluck_examples(glmdata,show_region)

[~,t0] = min(abs(glmdata.t_mids-0));
[~,t_end] = min(abs(glmdata.t_mids - 500));

% find units with 1 dominant predictor
ifsig = glmdata.sigunit(:,t0:t_end,3:4) & strcmp(glmdata.region,show_region);

idx_both = find(sum(sum(ifsig,2)>4,3)==2);
idx_dir = find(sum(ifsig(:,:,1),2)>7 & ...
    sum(ifsig(:,:,2),2)<4);
idx_val = find(sum(ifsig(:,:,1),2)<4 & ...
    sum(ifsig(:,:,2),2)>7);

% check CPD for these units
cpd = glmdata.cpd(:,t0:t_end,3:4);
keep_both = find(sum(max(cpd(idx_both,:,:),[],2)>.12,3)==2);
keep_dir = find(max(cpd(idx_dir,:,1),[],2)>.12);
keep_val = find(max(cpd(idx_val,:,2),[],2)>.12);

% example_units = [idx_both(keep_both);idx_dir(keep_dir);idx_val(keep_val)];

example_units = [idx_val(keep_val)];
end

