function show_cpd_heatmaps(glmdata, out_prefix)

region_names = unique(glmdata.region);
subject_names = unique(glmdata.subject);

nreg = length(region_names);
nsubj = length(subject_names);

useful_vars = glmdata.glmvars(2:4);

nmids = length(glmdata.t_mids);

% define colormap
X = [1 1 1;1 1 1 ; 1 1 1;.5 1 1; 0 .5 1;0 0 1; 0 0 .5;0 0 0;...
    .5 0 0; 1 0 0 ; 1 .5 0 ; 1 1 .5;1 1 1;1 1 1; 1 1 1 ;...
    1 1 1;1 1 1;.8 .8 .8; .6 .6 .6; .4 .4 .4; .2 .2 .2; 0 0 0];

X_smooth = mat2colormap(X);

if ~isempty(strfind(out_prefix,'alltrials'))
    fancy_label = {'','','';...
        'free','forced','not sig.';...
        'contra','ipsi','not sig.';...
        '+ maxval','- maxval','not sig.';...
        '','',''};
elseif ~isempty(strfind(out_prefix,'free'))
    fancy_label = {'','','';...
        'contra','ipsi','not sig.';...
        '+ maxval','- maxval','not sig.';...
        '+ \Deltaval','- \Deltaval','not sig.';...
        '','',''};
end

[~,t0] = min(abs(glmdata.t_mids-0));
[~,t_start] = min(abs(glmdata.t_mids - (-400)));
[~,t_end] = min(abs(glmdata.t_mids - 500));
[~,t_disp_end] = min(abs(glmdata.t_mids - 600));

for R = 1:nreg
    for S = 1:nsubj
        
        % restrict units
        region = region_names{R};
        subject = subject_names{S};
        
        idx = strcmp(glmdata.region,region_names{R}) & ...
            strcmp(glmdata.subject,subject_names{S});
        
        nunits = sum(idx);
        
        figure('Name',[subject,': ',region]);
        
        for V = 2:4
            
            % CPD for variable
            cpd = glmdata.cpd(idx,:,V);
            
            % get 1st sig bin t_pics>0 & beta sign
            
            ifsig = glmdata.sigunit(idx,t0:t_end,V);
            beta = glmdata.beta(idx,t0:t_end,V);
            firstbin = nan(nunits,1);
            betasign = nan(nunits,1);
            
            for u = 1:nunits
                frst = find(ifsig(u,:),1);
                if isempty(frst)
                    firstbin(u) = t0;
                    betasign(u) = 0;
                else
                    firstbin(u) = frst;
                    betasign(u) = 2*(beta(u,frst)>0)-1;
                end
            end
            
            % split unit sby beta sign
            pos_units = find(betasign>0);
            neg_units = find(betasign<0);
            zero_units = find(betasign==0);
            
            unit_counts = [length(pos_units),length(neg_units),length(zero_units)];
            unit_breaks = cumsum([0,unit_counts(1:2)]);
            
            [~,pos_idx] = sort(firstbin(pos_units),'descend');
            [~,neg_idx] = sort(firstbin(neg_units));
            [~,zero_idx] = sort(firstbin(zero_units));
            
            % get log10(cpd)
            cpd10 = log10(glmdata.cpd(idx,:,V));
            
            % recolor (to fit 1 colormap)
            rng = [-3 -1];
            cpd10(cpd10<rng(1)) = rng(1);
            cpd10(cpd10>rng(2)) = rng(2);
            
            cpd10_pos = -cpd10(pos_units(pos_idx),:)-1;
            cpd10_neg = cpd10(neg_units(neg_idx),:)+1;
            cpd10_zero = cpd10(zero_units(zero_idx),:)+5;
            
            % show!
            subplot(1,4,6-V);
            imagesc([cpd10_pos;cpd10_neg;cpd10_zero])
            caxis([-2 4])
            
            hold on
            
            % add t=0 line
            plot([t0 t0],.5+[0 nunits],'k','LineWidth',2)
            
            xlim([t_start,t_disp_end])
            
            % display correct x,y ticks
            temp = t0:-10:1;
            ticks = [temp(end:-1:2),t0:10:nmids];
            set(gca,'XTick',ticks,'XTickLabel',glmdata.t_mids(ticks)+1,...
                'YTick',0:100:nunits);
            
            % add unit count text & breaks
            for i = 1:3
                
                if i>0
                    plot(0.5+[0,nmids],unit_breaks(i)+[0.5 0.5],'k')
                end
                text(t0-13,unit_breaks(i) + 0.02*nunits,fancy_label{V,i})
                text(t0-13,unit_breaks(i) + 0.04*nunits,['n=',num2str(unit_counts(i))])
            end
            
            title(glmdata.glmvars{V})
            
            % axis labels
            xlabel('time from pics (ms)')
            ylabel('units')
        end
        
        subplot(1,4,1)
        imagesc([linspace(0,2,100)',linspace(0,-2,100)',linspace(4,2,100)'])
        hold on
        plot([1.5 2.5; 1.5 2.5],[.5 .5; 100.5 100.5],'k')
        set(gca,'YTick',[1,50,100],'YTickLabel',{'10%','1%','0.1%'},...
            'XTick',1:3,'XTickLabel',{'\beta>0','\beta<0','not sig.'},...
            'XTickLabelRotation',90,'XAxisLocation','top')
        caxis([-2 4])
        ylim([1 100])
        
        colormap(X_smooth)
        
        set(gcf,'Position',[100 100 900 800])
        
        % save fig
        print([out_prefix,'_',subject,'_',region,'.svg'],'-dsvg')
        
    end
end







end

