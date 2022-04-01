function [show_svgs, slice_images, corners_x, corners_z] = ...
    show_atlas(fname_atlas,region,AP_slices)

% load volume
nii = load_nii(fname_atlas);

vol = nii.img; % ML, AP, DV
dims = nii.hdr.dime.pixdim(2:4);
origin = nii.hdr.hist.originator(1:3);

% adjust (0,0,0) to earbars/eyebars
zero = round(origin-[0.5,22,18.5]./dims);

% cut down vol to useful display area
if strcmp(region,'CN')
    ML_range = [-100:-1]; % vox
    AP_range = [5:230]; % vox
    DV_range = [80:230]; % vox
    
elseif strcmp(region,'OFC')
    ML_range = [-100:-1]; % vox
    AP_range = [100:300]; % vox
    DV_range = [80:230]; % vox
end

vol_cube = vol(zero(1) + ML_range,...
    zero(2) + AP_range,...
    zero(3) + DV_range); % coords for atals ML, AP, DV

% which slices to show?
if isempty(AP_slices) % even spacing through AP_range (good when first starting, orienting)
    show_slices = 1:20:length(AP_range);
else
    show_slices = nan(size(AP_slices)); % show slices at specific APs (good for final figs)
    for j = 1:length(AP_slices)
        [~,idx] = min(abs(AP_range*dims(2) - AP_slices(j)));
        show_slices(j) = idx;
    end
end

% show slices!
corners_x = [min(ML_range), min(ML_range); max(ML_range), max(ML_range)]*dims(1);
corners_y = [0 0; 0 0];
corners_z = [min(DV_range), max(DV_range); min(DV_range), max(DV_range)]*dims(3);

slice_images = {};
for i = show_slices
    
    if strcmp(region,'CN')
        cnslice = permute(vol_cube(:,i,:)==82,[1,3,2]);
        putslice = permute(vol_cube(:,i,:)==98,[1,3,2]);
        
        slice = cnslice + putslice*2;
        
    elseif strcmp(region,'OFC')
        slice = double(permute(vol_cube(:,i,:)==216,[1,3,2]));
    end
    
    slice_images{end+1} = slice;
    
    x = surf(corners_x,corners_y+AP_range(i)*dims(2),corners_z,...
        'CData',slice,'FaceColor','texturemap','FaceAlpha',.3,...
        'BackFaceLighting','lit','AmbientStrength',1);
    
    % add other hemi! (optional)
    x = surf(-corners_x,corners_y+AP_range(i)*dims(2),corners_z,...
        'CData',slice,'FaceColor','texturemap','FaceAlpha',.3,...
        'BackFaceLighting','lit','AmbientStrength',1);
    
end

% find *.svg corresponding to slices
show_svgs = AP_range(show_slices)+zero(2);

end
