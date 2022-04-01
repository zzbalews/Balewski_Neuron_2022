function [AP_slices,track_coords_disc,track_coords_size] = add_units(ses_info,dir_hdd,dir_info,show_region)

subject_names = unique(ses_info{:,2});

% possible region names
ofc = {'ofc','OFC','LOFC','ROFC'};
cn = {'cn','CN','caudate','LCN','RCN'};
acc = {'acc','ACC','LACC','RACC'};

nses = size(ses_info,1);

coords = struct();
for S = 1:length(subject_names)
    coords.(subject_names{S}) = nan(0,3);
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
    
    % update by subject
    coords.(subject) = cat(1, coords.(subject), coords_ML_AP_DV(idx,:));
    
    
end

track_AP = [];
track_coords_disc = cell(1,2);
track_coords_size = cell(1,2);

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
    
    % disc to make nice plots
    if strcmp(subject,'Chap')
        ML_disc = make_disc(MLAPDV(:,1),[-11.5:11.5]);
    elseif strcmp(subject,'George')
        ML_disc = make_disc(MLAPDV(:,1),[-11:11]);
    end
    
    if strcmp(show_region,'CN')
        AP_disc = make_disc(MLAPDV(:,2),linspace(22,32,3));
    elseif strcmp(show_region,'OFC')
        AP_disc = make_disc(MLAPDV(:,2),[29,33,35,39]);%linspace(26,40,9));
    end
    track_AP = cat(2,track_AP,unique(AP_disc));
    
     if strcmp(show_region,'CN')
         DV_disc = make_disc(MLAPDV(:,3),[12:28]);
    elseif strcmp(show_region,'OFC')
         DV_disc = make_disc(MLAPDV(:,3),[12:1:28]);
    end
   
    
    % collapse by disc loc
    [coords_disc,~,idx] = unique([ML_disc',AP_disc',DV_disc'],'rows');
    track_coords_disc{S} = coords_disc;
    
    npts = size(coords_disc,1);
    
    coords_size = nan(npts,1);
    for n = 1:npts
        coords_size(n) = sum(idx==n);
    end
    track_coords_size{S} = coords_size;
    
    % plot
%     scatter3(MLAPDV(:,1),MLAPDV(:,2),MLAPDV(:,3),'filled','MarkerFaceAlpha',.3); % use raw coords
    scatter3(coords_disc(:,1),coords_disc(:,2),coords_disc(:,3),coords_size,'filled','MarkerFaceAlpha',.8); % use disc coords
end

AP_slices = unique(track_AP);

end

