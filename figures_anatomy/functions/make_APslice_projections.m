function make_APslice_projections(AP_slices, coords_disc, coords_size,...
    show_svgs, slice_images, corners_x, corners_z, dir_out,show_region)


for i = 1:length(AP_slices)
    
    figure; hold on
    daspect([1 1 1])
    
    % show anatomy
    x = surf(corners_x,corners_z,[0 0;0 0],...
        'CData',slice_images{i},'FaceColor','texturemap','FaceAlpha',.3,...
        'BackFaceLighting','lit','AmbientStrength',1);
    x = surf(-corners_x,corners_z,[0 0;0 0],...
        'CData',slice_images{i},'FaceColor','texturemap','FaceAlpha',.3,...
        'BackFaceLighting','lit','AmbientStrength',1);
    
    if strcmp(show_region,'CN')
        scale = 0.5;
    elseif strcmp(show_region,'OFC')
        scale = 2;
    end
    
    % show units
    for s = 1:2
        idx = coords_disc{s}(:,2)==AP_slices(i);
        scatter(coords_disc{s}(idx,1),coords_disc{s}(idx,3),5+coords_size{s}(idx)*scale,'filled')
    end
    
    % add size refs
     if strcmp(show_region,'CN')
          scatter([10,10],[29,28],5+[10,100]*scale,'filled','MarkerFaceColor',[0 0 0])
    text([11,11],[29,28],{'10','100'})
    elseif strcmp(show_region,'OFC')
          scatter([10,10],[29,28],5+[5,50]*scale,'filled','MarkerFaceColor',[0 0 0])
    text([11,11],[29,28],{'5','50'})
     end
    
  
    
    % format
    if strcmp(show_region,'CN')
        xlim([-15 15]);
        ylim([12 32]);
    elseif strcmp(show_region,'OFC')
        xlim([-15 15]);
        ylim([12 32]);
    end
    
    set(gcf,'Position',[50 50 700 500])
    

    title({show_region,['svg #',num2str(show_svgs(i))],['AP ',num2str(AP_slices(i))]})
end

end

