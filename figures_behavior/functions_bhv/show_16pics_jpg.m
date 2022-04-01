function show_16pics_jpg(bhvdata,dir_stims,which_img)
%show jpgs in 4x4 pic arrangement

%% show subject results (separte plots)
figure;

ranking = {'L','ML','MH','H'};
px = 100;

subject_names = unique(bhvdata.subject);
for s = 1:length(subject_names)
    
    pic_matrix = zeros(4*px,4*px,3);
    
    % restrict to ('trialtype') trials from subject
    subject = subject_names{s};
    tr_subj = strcmp(bhvdata.subject,subject);
    
    % get amnt/prob levels
    amnts = unique(bhvdata.amnt(tr_subj,:));
    amnts = amnts(amnts>0);
    probs = unique(bhvdata.prob(tr_subj,:));
    probs = probs(probs>0);
    
    % load all images, arrange correctly
    for a = 1:4
        for p = 1:4
            
            if strcmp(which_img,'jpg')
            fname = ['amnt',ranking{a},'_prob',ranking{p},'.jpg'];
            temp = imread(fullfile(dir_stims,[subject,'_stims'],fname));
            
            else
               idx = find(bhvdata.amnt(:,1)==amnts(a) & ...
                   bhvdata.prob(:,1)==probs(p) & ...
                   tr_subj,1);
               breaks = [30,80,160,240];
               temp = breaks(bhvdata.valbin_expval(idx,1))*ones(px,px,3);
                
            end
            idx_row = px*a + (0:-1:(1-px));
            idx_cols = px*p + ((1-px):0);
            
            pic_matrix(idx_row,idx_cols,:) = temp;
            
        end
    end
    
    % show imagesc
    subplot(1,2,s)
    
    imagesc(pic_matrix/256);
    
    hold on
    daspect([1 1 1]);
    
    for i = 100.5:100:300.5
        plot([-3,404],[i i],'k','LineWidth',2)
        plot([i i],[-3,404],'k','LineWidth',2)
    end
    
    title(subject)
    set(gca,'YDir','normal',...
        'XTick',50:100:400,'XTickLabel',round(probs,2),...
        'YTick',50:100:400,'YTickLabel',round(amnts,1))
    xlabel('reward probability')
    ylabel('juice amount')
   
end


end

