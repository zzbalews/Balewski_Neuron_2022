function show_unit_profiles(glmfile,which_trials,t_range)

Nfiles = length(glmfile);




if strcmp(which_trials,'alltrials')
clrs = [238,204,136;240,144,89;184,63,95]/256;
labels = {'trial type','lever dir','max subj val'};
keep_vars = 2:4;
elseif strcmp(which_trials,'free')
    clrs = [240,144,89;184,63,95;173, 85, 170]/256;
    labels = {'lever dir','max subj val','\Delta subj val'};
    keep_vars = 2:4;
elseif strcmp(which_trials,'freesimple')
    labels = {'lever dir','max subj val'};
    clrs = [240,144,89;184,63,95]/256;
    keep_vars = 2:3;
else
    error('trial type invalid')
end

nV = length(keep_vars);
%% load glm output & session details

unit_names = {};
sigunit = [];
unit_fr = [];
for f = 1:Nfiles

    temp = load(glmfile{f});
    
    unit_names = cat(1,unit_names,temp.output.unit_names);
    ses_name = strrep(temp.output.session,'_',' ');
    
    if Nfiles==1
        region = temp.meta.region;
    else
        region = temp.meta.region(2:end);
    end
    

    if f==1
        varnames = temp.output.varnames(keep_vars);
        t_mids = temp.t_mids;
    end
    
    sigunit = cat(1,sigunit,temp.output.sigunit);
    unit_fr = cat(1,unit_fr,temp.output.unit_fr');
end

nunits = length(unit_names);
%% show figure
figure; 

% # sig units 0-500ms
subplot(2,1,1);hold on

idx = t_mids > t_range(1) & t_mids < t_range(2);
nsig = permute(sum(sum(sigunit(:,idx,keep_vars),2)>0),[1,3,2]);

yrng = [0,1.2*max(nsig)];
for i = 1:nV
    bar(i,nsig(i),'FaceColor',clrs(i,:))
    text(i,nsig(i) + yrng(2)*0.05,[num2str(round(100*nsig(i)/nunits)),'%'],'HorizontalAlignment','center');
end

set(gca,'XTick',1:nV,...
    'XTickLabel',labels,...
    'XTickLabelRotation',-90)

ylabel(['# sig units (',num2str(t_range(1)),' : ',num2str(t_range(2)),' ms)'])
ylim(yrng)

title({ses_name,[region,': ',which_trials,' (',num2str(nunits),' units)']})

% firing rates
subplot(2,1,2); hold on
histogram(unit_fr,'BinWidth',1,'FaceColor',[.5 .5 .5])
xlabel('firing rate (Hz)')
ylabel('count')


% format & save
set(gcf,'Position',[100 100 250 550])

temp = strsplit(glmfile{1},'_');
temp{end} = [region,'.png'];
figname = strjoin(temp,'_');

print(figname,'-dpng')

end

