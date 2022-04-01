function show_psychometric(bhvdata)
%plot psychometric fits for each subject

% colors
clr_saccL = [80, 180, 66]/256;
clr_saccR = [160, 94, 179]/256;

%% show subject results (separte plots)
figure;
subject_names = unique(bhvdata.subject);
for s = 1:length(subject_names)
    
    % restrict to ('trialtype') trials from subject
    subject = subject_names{s};
    tr_subj = strcmp(bhvdata.subject,subject) & ...
        bhvdata.trialtype==2 & ...
        bhvdata.lever~=0 & ~isnan(bhvdata.saccloc(:,1));
    
    
    % do psycho fit for model
    lever = bhvdata.lever(tr_subj);
    jpg = bhvdata.jpg(tr_subj,:);
    amnt = bhvdata.amnt(tr_subj,:);
    prob = bhvdata.prob(tr_subj,:);
    firstsacc = bhvdata.saccloc(tr_subj,1);
    
    [w_names,w_fit,BIC,R2] = ...
        fit_discount_sacc('expval',lever,jpg,amnt,prob,firstsacc);
    
    vdiff = -6:.1:6;
    goright_fit_saccL = ...
        1-softmax_discount_sacc('expval',{vdiff},-1,w_fit);
    goright_fit_saccR = ...
        1-softmax_discount_sacc('expval',{vdiff},1,w_fit);
    
    % plot model
    subplot(1,2,s); hold on
    plot([0 0],[0 1],'k:')
    plot([-6 6],[.5 .5],'k:')
    
   
    modL = plot(-vdiff,goright_fit_saccL,'LineWidth',2,'Color',[clr_saccL]);
    modR = plot(-vdiff,goright_fit_saccR,'LineWidth',2,'Color',[clr_saccR]);
    
    % plot choice data
    tr_vdiff = diff(bhvdata.subjval_expval(tr_subj,:),[],2);
    maxdiff = ceil(max(abs(tr_vdiff)));
    
    edges = linspace(-maxdiff,maxdiff,25);
    midpoints = mean([edges(1:end-1); edges(2:end)]);
    
    tr_vdiff_disc = midpoints(discretize(tr_vdiff,edges))';
    
    goright_data_saccL = nan(length(midpoints),2);
    goright_data_saccR = nan(length(midpoints),2);
    
    for m = 1:length(midpoints)
        % sacc1 = left
        idx = tr_vdiff_disc==midpoints(m) & firstsacc==-1;
        if sum(idx)==0
            goright_data_saccL(m,:) = [NaN,NaN];
        else
            goright_data_saccL(m,:) = [mean(lever(idx)==1), sum(idx)];
        end
        % sacc1 = right
        idx = tr_vdiff_disc==midpoints(m) & firstsacc==1;
        if sum(idx)==0
            goright_data_saccR(m,:) = [NaN,NaN];
        else
            goright_data_saccR(m,:) = [mean(lever(idx)==1), sum(idx)];
        end
    end
    
    datL = scatter(midpoints(1:10),goright_data_saccL(1:10,1),10+goright_data_saccL(1:10,2)/10,...
        'filled','MarkerFaceColor',clr_saccL,'MarkerFaceAlpha',.7);
    
    datR = scatter(midpoints,goright_data_saccR(:,1),10+goright_data_saccR(:,2)/10,...
        'filled','MarkerFaceColor',clr_saccR,'MarkerFaceAlpha',.7);
    
    scatter(midpoints(11:end),goright_data_saccL(11:end,1),10+goright_data_saccL(11:end,2)/10,...
        'filled','MarkerFaceColor',clr_saccL,'MarkerFaceAlpha',.7);
    
    scatter([4 4 4],[.4 .3 .2],5+[50 200 1000]/10,'filled','MarkerFaceColor',[.5 .5 .5])
    text([4.3 4.3 4.3],[.4 .3 .2],{'50','200','1000'},'HorizontalAlignment','left')
    
    title(subject)
    xlabel('value_{right} - value_{right}')
    ylabel('proportion choose right')
    
    legend([datL,datR,modL,modR],{'1st sacc Left','1st sacc Right','model','model'},...
        'box','off','Location','northwest')
    
    
    %% summary stats
    disp(subject)
    disp(['... model exp var R2 = ',num2str(round(100*R2,2)),'%'])
    
    correct_lever = double(tr_vdiff>0);
    correct_lever(correct_lever==0) = -1;
    disp(['... % trials chose high>low = ',num2str(mean(lever==correct_lever))])
    
    disp(['... % trials 1st sacc to high = ',num2str(mean(firstsacc==correct_lever))])
    
end

end

