function add_psychometric(cleanfile,str_model_name)
%add_psychometric(cleanfile,str_model_name)
%   add subjective value & value bin of pictures
%   cleanfile = processed bhv file
%   str_model_name = model to use for value
%   ses_fullname = string of session ID for figures

%% load processed behavior
load(cleanfile)

% skip if subj value alread exists
if ismember('subjval',fieldnames(trialinfo))
    return
end
%% get trial info

% use only completed free trials
keep = trialinfo.trialtype==2 & ~isnan(trialinfo.lever) & ~isnan(trialinfo.firstsaccloc);
    
% choice direction
lever = trialinfo.lever(keep); % -1=left, +1=right

% trial info (column 1 = left, 2 = right)
jpg = trialinfo.jpg(keep,:); % picture ID (17=empty space)
amnt = trialinfo.amnt(keep,:); % juice amount, in pump Voltage
prob = trialinfo.prob(keep,:); % prob of juice reward
firstsacc = trialinfo.firstsaccloc(keep); % -1=left, +1=right

%% do fit, according to model
if strcmp(str_model_name,'linear discount + first saccade')
    [w_names,w_fit,BIC,R2] = fit_discount_sacc('lin',lever,jpg,amnt,prob,firstsacc);
    subjval = lindiscount(trialinfo.amnt,trialinfo.prob,w_fit(1));
else
    disp('ERROR: psychometric model name not recognized')
    return
end

%% update processed behavior
temp = trialinfo.jpg==17; % empty spaces on forced
subjval(temp) = NaN;
trialinfo.subjval = subjval;
trialinfo.valbin = bin_val(subjval,4);

save(cleanfile,'trialinfo')


end

