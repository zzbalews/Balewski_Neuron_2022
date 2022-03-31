function probeID = get_contact_location(subject,...
    probeID,recinfo,transform)

space = 0.75; % mm distance between grid holes

switch subject
    case 'Chap'
        gridwidth = 12.5; % mm thickness of recording grid
        turn_dist = 1/3; % mm traveled per screw turn
    case 'George'
        gridwidth = 14; % mm thickness of recording grid
        turn_dist = 1/4; % mm traveled per screw turn
end

% load probe coord transforms
load(transform)

% restrict to successful probes
goodIDs = unique(probeID{:,3});

figure; hold on
daspect([1 1 1])
xlabel('ML')
ylabel('AP')
zlabel('DV')

% track coords
for g = 1:length(goodIDs)
    
    % probe #
    prnum = goodIDs(g);
    idx = find(recinfo{:,1}==prnum);
    
    % probe details
    coords = recinfo{idx,3:4};
    hemi = recinfo{idx,2}{1};
    region = recinfo{idx,5}{1};
    
    temp = strsplit(recinfo{idx,8}{1},'-');
    chan = str2num(temp{1}):str2num(temp{2}); % channels on this probe
    
    contact_dist = recinfo{idx,14}/1000; % mm spacing between contacts
    
    guide_mm = recinfo{idx,9};
    turns = recinfo{idx,10};
    
    tiltdeg = recinfo{idx,11:12}; % ML, AP in deg
    
    % adjust trajectories by ML & AP tilts
    ML_elec_rad = deg2rad(tiltdeg(1));
    tilt_ML = [
        cos(ML_elec_rad) 0 sin(ML_elec_rad);...
        0 1 0;...
        -sin(ML_elec_rad) 0 cos(ML_elec_rad)];
    
    AP_elec_rad = deg2rad(tiltdeg(2));
    tilt_AP = [
        1 0 0;...
        0 cos(-AP_elec_rad) -sin(-AP_elec_rad);...
        0 sin(-AP_elec_rad) cos(-AP_elec_rad)];
    
    % which region?
    if strcmp(region,'ofc') | strcmp(region,'OFC')
        clr = [0.7 0 0.7];
    else
        clr = [0 0.7 0];
    end
    
    % which hemi transform?
    switch hemi
        case 'left'
            move = left;
        case 'right'
            move = right;
    end
    
    % starting grid hole
    switch subject
        case 'Chap'
            hole2D =  [coords(1),coords(2),0]*space;
        case 'George'
            hole2D = [-coords(1),-coords(2),0]*space;
    end
    
    hole3D = hole2D * move.rot_1 * move.rot_2 * move.rot_3 +...
        move.shift;
    
    scatter3(hole3D(1),hole3D(2),hole3D(3),'filled','k')
    
    % probe end
    extended_mm = gridwidth + guide_mm + turns*turn_dist;
    stop2D = hole2D + [0 0 -extended_mm];
    
    stop3D = stop2D * move.rot_1 * move.rot_2 * move.rot_3 * ...
        tilt_ML * tilt_AP + ...
        move.shift;
    
    temp = [stop3D; hole3D];
    plot3(temp(:,1),temp(:,2),temp(:,3),'Color',clr)
    
    % add channel holes
    switch subject
        case 'Chap'
            temp = stop3D-hole3D;
        case 'George'
            temp = hole3D-stop3D;
    end
    unitnorm = temp./sqrt(sum(temp.^2,2));
    
    contacts = -(0.3 + contact_dist*(max(chan)-chan)); % dist from tip of probe
    markers3D = repmat(stop3D,length(contacts),1) + contacts'*unitnorm;
    
    scatter3(markers3D(:,1),markers3D(:,2),markers3D(:,3),5,'filled','k')
    
    % update with coords
    for i = 1:length(chan)
        idx = find(probeID{:,1}==chan(i));
        probeID(idx,'ML') = {markers3D(i,1)};
        probeID(idx,'AP') = {markers3D(i,2)};
        probeID(idx,'DV') = {markers3D(i,3)};
    end
end
view([150 25])

end
