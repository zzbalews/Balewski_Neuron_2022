function show_16pics_heatmap(bhvdata,variable,trialtype)
%show 'variable' in 4x4 pic arrangement

% restrict to selected trial type
switch trialtype
    case 'free'
        tr_trtype = bhvdata.trialtype==2 & ...
            bhvdata.lever~=0;
    case 'forced'
        tr_trtype = bhvdata.trialtype==1 & ...
            bhvdata.lever~=0;
end

% map lever options: -1 = left, +1 = right
lever_opts = [-1,1];

%% show subject results (separte plots)
figure('Name',[trialtype,': ',variable]);

subject_names = unique(bhvdata.subject);
for s = 1:length(subject_names)
    
    % restrict to ('trialtype') trials from subject
    subject = subject_names{s};
    tr_subj = strcmp(bhvdata.subject,subject) & tr_trtype;
    
    % get amnt/prob levels
    amnts = unique(bhvdata.amnt(tr_subj,:));
    amnts = amnts(amnts>0);
    probs = unique(bhvdata.prob(tr_subj,:));
    probs = probs(probs>0);
    
    % get all trials with each picture
    track_ntr = zeros(4);
    track_variable = zeros(4);
    
    for a = 1:length(amnts)
        for p = 1:length(probs)
            for idx = 1:2 % left/right
                
                subset = tr_subj & ...
                    bhvdata.amnt(:,idx)==amnts(a) & ...
                    bhvdata.prob(:,idx)==probs(p);
                
                track_ntr(a,p) = track_ntr(a,p) + ...
                    sum(subset);
                
                switch variable
                    case 'choice'
                        % show prop chosen
                        track_variable(a,p) = track_variable(a,p) + ...
                            sum(bhvdata.lever(subset)==lever_opts(idx));
                        
                    case 'RT'
                        % show lever response time
                        track_variable(a,p) = track_variable(a,p) + ...
                            sum(bhvdata.rt(subset));
                        
                    case 'saccloc'
                        % show 1st sacc location
                        track_variable(a,p) = track_variable(a,p) + ...
                            sum(bhvdata.saccloc(subset,1)==lever_opts(idx));
                        
                    case 'sacctime'
                        % show 1st sacc time
                        track_variable(a,p) = track_variable(a,p) + ...
                            nansum(bhvdata.MLsacc(subset,1));
                        
                        track_ntr(a,p) = track_ntr(a,p) - ...
                            sum(isnan(bhvdata.MLsacc(subset,1)));
                        
                end
            end
        end
    end
    
    subplot(1,2,s)
    
    imagesc(track_variable ./ track_ntr);
    daspect([1 1 1]);
    colorbar
    title(subject)
    set(gca,'YDir','normal',...
        'XTick',1:4,'XTickLabel',round(probs,2),...
        'YTick',1:4,'YTickLabel',round(amnts,1))
    xlabel('reward probability')
    ylabel('juice amount')
    
    switch variable
        case 'choice'
            caxis([0 1])
        case 'saccloc'
            caxis([0 1])
        case 'sacctime'
            caxis([150 375])
    end
    
end

colormap(hot)
end

