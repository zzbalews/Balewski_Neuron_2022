function glmdata = load_glmresults(ses_info, dir_data, dir_unitglms, which_model)

% variables
varnames = {'experiment','session','subject','unit_names','region','hemisphere',...
    'beta','cpd','sigunit'};

% init glmdata struct
glmdata = struct();
for v = 1:length(varnames)
   glmdata.(varnames{v}) = []; 
end

nses = size(ses_info,1);

vifs = nan(nses,1);

for s = 1:nses
   
    % get session details
    experiment = ses_info{s,1}{1};
    subject = ses_info{s,2}{1};
    num = ses_info{s,3}{1};
    date = strrep(ses_info{s,4}{1},'/','');
    
    session = strjoin({num,date},'_');
    
    % get session glm output name
    fileloc = fullfile(dir_data,experiment,subject,dir_unitglms,session);
    temp = dir([fileloc,'/*',which_model,'*.mat']);
    filenames = {temp.name};
    
    % read in!
    for f = 1:length(filenames)
        
        % load file
        load(fullfile(fileloc,filenames{f}),'meta','output','t_mids');
        
        % get hemi & region
        hemi = meta.region(1);
        reg = meta.region(2:end);
        
        % save important vars
        if f==1
            vifs(s) = output.vif;
        end
        
        nunits = length(output.unit_names);
        
        output.experiment = repmat({experiment},nunits,1);
        output.subject = repmat({subject},nunits,1);
        output.session = repmat({session},nunits,1);
        output.region = repmat({reg},nunits,1);
        output.hemisphere = repmat({hemi},nunits,1);
        
        for v = 1:length(varnames)
           glmdata.(varnames{v}) = cat(1, glmdata.(varnames{v}), output.(varnames{v})); 
        end       
       
    end
end

glmdata.t_mids = t_mids;
glmdata.glmvars = output.varnames;

end

