function make_nice_psths(bhvdata,glmdata,example_units,dir_data,dir_pl2)

using_sessions = unique(glmdata.session(example_units));

t_mids = glmdata.t_mids;


figure;
for ses = 1:length(using_sessions)
    
    tic
    % get bhv for this session
    tr = strcmp(bhvdata.session,using_sessions{ses});
    mini = struct();
    mini.maxval = max(bhvdata.valbin_expval(tr,:),[],2);
    mini.lever = bhvdata.lever(tr);
    
   % check out all units in this session
   unit_counter = find(strcmp(glmdata.session(example_units),using_sessions{ses}))
   current_units = example_units(unit_counter);
   
   % get spk filenames
   experiment = glmdata.experiment{current_units(1)};
   subject = glmdata.subject{current_units(1)};
    session = glmdata.session{current_units(1)};
    current_units_names = {glmdata.unit_names{current_units}};
   current_units_region = {glmdata.region{current_units}};
   current_units_hemi = {glmdata.hemisphere{current_units}};
   
   
    fileloc = fullfile(dir_data,experiment,subject,dir_pl2,session);
    temp = dir([fileloc,'/*w100_s25_pics*.mat']);
    fnames = {temp.name};
    
 
    % check all files, load corresponding spk row
    current_units_fr = [];
    for f = 1:length(fnames)
        
    load(fullfile(fileloc,fnames{f}),'SPKfr_names');
    
    for u = 1:length(current_units)
       
        rowid = find(strcmp(SPKfr_names,current_units_names{u}));
        if isempty(rowid)
            continue
        else
            ff = matfile(fullfile(fileloc,fnames{f}));
            current_units_fr(:,:,u) = ff.SPKfr_units(:,:,rowid);
        end
    end
        
    end
    
    % if sig?
    current_units_ifsigdir = glmdata.sigunit(current_units,:,3);
    current_units_ifsigmaxval = glmdata.sigunit(current_units,:,4);

    % show all unit psths
    
    for u = 1:length(current_units)
                
        % plot contra/ipis
        ax1 = subplot(2,5,unit_counter); hold on
        
        if strcmp(current_units_hemi{u},'L')
            tr_contra = mini.lever==1;
            tr_ipsi = mini.lever==-1;
        elseif strcmp(current_units_hemi{u},'R')
            tr_contra = mini.lever==-1;
            tr_ipsi = mini.lever==1;
        end
        
        
%         plot(t_mids,mean(current_units_fr(tr_ipsi,:,u)),...
%             'Color',.2*[255,215,0]/256,'LineWidth',2);
%         plot(t_mids,mean(current_units_fr(tr_contra,:,u)),...
%             'Color',.9*[255,215,0]/256,'LineWidth',2);
%         
        shadedErrorBar(t_mids, mean(current_units_fr(tr_ipsi,:,u)),...
            std(current_units_fr(tr_ipsi,:,u))./sqrt(sum(tr_ipsi)),...
            'lineprops',{'Color',[0 0 0],'LineWidth',1})
        
        shadedErrorBar(t_mids, mean(current_units_fr(tr_contra,:,u)),...
            std(current_units_fr(tr_contra,:,u))./sqrt(sum(tr_contra)),...
            'lineprops',{'Color',[255,215,0]/256,'LineWidth',1})
        
        
        y = get(gca,'YLim');
        temp = y(2)-0.1*diff(y);
        
        ifsig = current_units_ifsigdir(u,:); 
        ifsig(ifsig==0) = NaN;
        plot(t_mids,ifsig.*zeros(size(t_mids)),'r','LineWidth',2)
        
        title({current_units_names{u},current_units_region{u},subject})
        
        % plot max value
        ax2 = subplot(2,5,unit_counter+5); hold on
        
        v1 = mini.maxval==1 & mini.lever~=0;
        v2 = mini.maxval==2 & mini.lever~=0;
        v3 = mini.maxval==3 & mini.lever~=0;
        v4 = mini.maxval==4 & mini.lever~=0;
        
        clrs = [linspace(0,.3,4)',linspace(0,.7,4)',linspace(0,1,4)'];
        
%         plot(t_mids,mean(current_units_fr(v1,:,u)),...
%             'Color',clrs(1,:),'LineWidth',2);
%         plot(t_mids,mean(current_units_fr(v2,:,u)),...
%             'Color',clrs(2,:),'LineWidth',2);
%         plot(t_mids,mean(current_units_fr(v3,:,u)),...
%             'Color',clrs(3,:),'LineWidth',2);
%         plot(t_mids,mean(current_units_fr(v4,:,u)),...
%             'Color',clrs(4,:),'LineWidth',2);
%         
        shadedErrorBar(t_mids,mean(current_units_fr(v1,:,u)),...
            std(current_units_fr(v1,:,u))./sqrt(sum(v1)),...
            'lineprops',{'Color',clrs(1,:),'LineWidth',1});
        
         shadedErrorBar(t_mids,mean(current_units_fr(v2,:,u)),...
            std(current_units_fr(v2,:,u))./sqrt(sum(v2)),...
            'lineprops',{'Color',clrs(2,:),'LineWidth',1});
        
        
         shadedErrorBar(t_mids,mean(current_units_fr(v3,:,u)),...
            std(current_units_fr(v3,:,u))./sqrt(sum(v3)),...
            'lineprops',{'Color',clrs(3,:),'LineWidth',1});
        
        
         shadedErrorBar(t_mids,mean(current_units_fr(v4,:,u)),...
            std(current_units_fr(v4,:,u))./sqrt(sum(v4)),...
            'lineprops',{'Color',clrs(4,:),'LineWidth',1});
        
        
        
        y = get(gca,'YLim');
        temp = y(2)-0.1*diff(y);
        
        ifsig = current_units_ifsigmaxval(u,:); 
        ifsig(ifsig==0) = NaN;
        plot(t_mids,ifsig.*zeros(size(t_mids)),'r','LineWidth',2)
        
                title({current_units_names{u},current_units_region{u},subject})

        % format
        linkaxes([ax1,ax2],'y')
        for i = [0,5]
        subplot(2,5,unit_counter+i)
        xlim([-300 600])
        y = get(gca,'YLim');
        plot([0 0],y+[0 5],'k')
        set(gca,'YTick',0:10:y(2),'XTick',[0,400])
        xlabel('time from pics (ms)')
        ylabel('firing rate (Hz)')
        end
       
        
        
    end
        
%        set(gcf,'Position',[100 100 185 400]) 
    
    toc
    end
    
    

end

