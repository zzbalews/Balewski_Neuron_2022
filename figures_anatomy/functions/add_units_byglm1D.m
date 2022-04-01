function [AP_slices,track_coords_disc,track_coords_size] = add_units_byglm1D(...
    ses_info,dir_hdd,dir_info,show_region,glmdata,keepdim)

subject_names = unique(ses_info{:,2});

% possible region names
ofc = {'ofc','OFC','LOFC','ROFC'};
cn = {'cn','CN','caudate','LCN','RCN'};
acc = {'acc','ACC','LACC','RACC'};

nses = size(ses_info,1);


%% load units in region, by subject

coords = struct();
unitID = struct();
unitGLM = struct();
for S = 1:length(subject_names)
    coords.(subject_names{S}) = nan(0,3);
    unitID.(subject_names{S}) = cell(0,2);
    unitGLM.(subject_names{S}).ifsig = nan(0,3);
end

for s = 1:nses
    
    % get session details
    experiment = ses_info{s,1}{1};
    subject = ses_info{s,2}{1};
    num = ses_info{s,3}{1};
    date = strrep(ses_info{s,4}{1},'/','');
    
    session = strjoin({num,date},'_');
    
    % load spk location file
    fileloc = fullfile(dir_hdd,subject,experiment,dir_info,session);
    temp = dir([fileloc,'/*_unit_MLAPDVloc.mat']);
    fname = temp(1).name;
    
    load(fullfile(fileloc,fname));
    
    
    if sum(coords_ML_AP_DV(:,2)<20)>0
        disp(session)
    end
    
    if sum(isnan(coords_ML_AP_DV(:,1)))>0
        disp(session)
    end
    
    
    % keep units in region
    if strcmp(show_region,'OFC')
        idx = ismember(region,ofc);
    elseif strcmp(show_region,'CN')
        idx = ismember(region,cn);
    elseif strcmp(show_region,'ACC')
        idx = ismember(region,acc);
    else
        error('do not recognize this region')
    end
    
    % update by subject:
    % unit coords
    coords.(subject) = cat(1, coords.(subject), coords_ML_AP_DV(idx,:));
    % unit id
    temp = [repmat({session},sum(idx),1),unit_names(idx)];
    unitID.(subject) = cat(1, unitID.(subject), temp);
    
end

%% get glm results

% determine if sig by term
t_mids = glmdata.t_mids;
t_range = t_mids>=0 & t_mids<=500;

temp = permute(sum(glmdata.sigunit(:,t_range,:),2)>0,[1,3,2]);
ifsig = [temp(:,3:4),temp(:,3)&temp(:,4)];
terms = [glmdata.glmvars(3:4),'both'];
terms_nice = {'lever direction','max value','lever * value'};
for S = 1:length(subject_names)
    subject = subject_names{S};
    nunits = length(coords.(subject));
    
    for u = 1:nunits
        
        glmidx = strcmp(glmdata.session, unitID.(subject){u,1}) & ...
            strcmp(glmdata.unit_names, unitID.(subject){u,2});
        
        unitGLM.(subject).ifsig(u,:) = ifsig(glmidx,:);
        
    end
end

%% to visualize, round positions to slices
track_AP = [];
track_coords_disc = cell(1,2);
track_coords_size = cell(1,2);

marker_shape = {'d','o'};

keepdim='DV';
figure; hold on
for S = 1:length(subject_names)
    subject = subject_names{S};
    MLAPDV = coords.(subject);
    
    % hacky fixes to align to standard brain...
    if strcmp(show_region,'CN')
        if strcmp(subject,'Chap')
            MLAPDV(:,1) = MLAPDV(:,1)-2;
            MLAPDV(:,2) = MLAPDV(:,2)-2;
            MLAPDV(:,3) = MLAPDV(:,3) + 1;
        elseif strcmp(subject,'George')
            idx = MLAPDV(:,1)<0;
            MLAPDV(idx,1) = MLAPDV(idx,1)-1;
            MLAPDV(~idx,1) = MLAPDV(~idx,1)+1;
            MLAPDV(:,2) = MLAPDV(:,2)+1;
            MLAPDV(:,3) = MLAPDV(:,3) + 7;
        end
        
    elseif strcmp(show_region,'OFC')
        if strcmp(subject,'Chap')
            0;
        elseif strcmp(subject,'George')
            MLAPDV(:,2) = MLAPDV(:,2) + 2;
            MLAPDV(:,3) = MLAPDV(:,3) + 6;
        end
    end
    
    MLAPDV(:,1) = abs(MLAPDV(:,1));
    
    % ignore 'collapse' dim
    if strcmp(keepdim, 'ML')
        dims_idx = 1;
%         dims_idx = [2,3];
%         dims_names = {'AP','DV'};
%         dims_lim = [21 33; 18 30];
    elseif strcmp(keepdim, 'AP')
        dims_idx = 2;
%         dims_idx = [1,3];
%         dims_names = {'ML','DV'};
%         dims_lim = [2 9; 18 30];
    elseif strcmp(keepdim,'DV')
        dims_idx = 3;
%         dims_idx = [1,2];
%         dims_names = {'ML','AP'};
%         dims_lim = [2 9; 21 34];
        
    end
    
    % disc to make nice plots, using flattened 2D projection
    flat_raw = MLAPDV(:,dims_idx);
    flat_disc = nan(size(flat_raw));
    for i = 1:size(flat_raw,2)
        low = floor(min(flat_raw(:,i)));
        high = ceil(max(flat_raw(:,i)));
        
        if S==2 & i==2
            low = low-0.5;
        end
        
        flat_disc(:,i) = make_disc(flat_raw(:,i), [low:.5:high]);
    end
    
    % group by disc location
    [coords_disc, ~, idx] = unique(flat_disc, 'rows');
    npts = size(coords_disc,1);
    
    coords_size = nan(npts,1);
    coords_sig = nan(npts,3);
    for n = 1:npts
        coords_size(n) = sum(idx==n);
        coords_sig(n,:) = mean(unitGLM.(subject).ifsig(idx==n,:),1);
    end
    
    %% do regression!
    Y = logit(coords_sig);
    Y_norm = zscore(Y);
    
    X = coords_disc;
    vif = max(diag(inv(corrcoef(X))));
    X_norm = [ones(npts,1), zscore(X)];
    
    [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X_norm, Y_norm);
    
    model_p
%     
%     
%     %% visualize
%     clrs_max = 0*[1 1 1; 1 1 1; 1 1 1]/256;
%     for j = 1:3
%         subplot(1,3,j); hold on
%         clrs = (1-coords_sig(:,j)) + coords_sig(:,j)*clrs_max(j,:);
%         
%         scatter(coords_disc(:,1),coords_disc(:,2),...
%             2*coords_size,clrs, 'filled',...
%             'MarkerFaceAlpha',1,'Marker',marker_shape{S},...
%             'MarkerEdgeAlpha',0);
%         
%         daspect([1 1 1])
%         
%         xlim(dims_lim(1,:))
%         ylim(dims_lim(2,:))
%         
%         view([0 -90])
%         title(terms_nice{j})
%         
%         text(dims_lim(1,1)+1,dims_lim(2,2)-(S/2),[subject(1),': model p=',num2str(round(model_p(j),3))])
%         
%         if model_p(j)<0.01
%             
%             disp([terms_nice{j},': ',subject])
%             disp(['   ',dims_names{1},' B=',num2str(round(B_hat(2,j),2)),', p=',num2str(round(coeff_p(2,j),3))])
%             disp(['   ',dims_names{2},' B=',num2str(round(B_hat(3,j),2)),', p=',num2str(round(coeff_p(3,j),3))])
%             
%         end
%         
%         xlabel(dims_names{1})
%         ylabel(dims_names{2})
%         
%     end
%     %      clrmap = repmat(1-[0:.01:1]',1,3);
%     %         colormap(clrmap)
%     %         colorbar
%     
%     
%     
    
end
% set(gcf,'Position',[200 200 1200 600])


end

