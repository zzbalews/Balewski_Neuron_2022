function show_avg_ppd_bylever(minidata,outdir,spksmooth)

%% load saved free post prob % delta
temp = strsplit(spksmooth,'_');
session = strjoin(temp(1:3),' ');
region = strrep(temp{end},'.mat','');

dataprefix = strjoin({outdir,[strrep(session,' ','_'),'_',region,'_value_free']},'/');
load([dataprefix,'.mat'])

if strcmp(region(1),'R')
    hemi = 'right';
    tt = {'contra','ipsi'};
elseif strcmp(region(1),'L')
    hemi = 'left';
    tt = {'ipsi','contra'};
end

%% split trials by left/right lever

L_tr = minidata.lever(free_tr,:)==-1;
R_tr = minidata.lever(free_tr,:)==1;

%% show post prob delta split by left/right trials
figure;
for i = 1:2
subplot(1,2,i); hold on
plot([0 0],[-50 200],'k')
plot([-2000,2000],[0,0],'k')
end

subplot(1,2,1)
a=plot(t_mids, mean(ch_ppd(L_tr,:)), 'Color', [1 0 0]);
b=plot(t_mids, mean(unch_ppd(L_tr,:)), 'Color', [0 0 1]);
legend([a,b],{'ch','unch'},'box','off','Location','northeast')

subplot(1,2,2)
a=plot(t_mids, mean(ch_ppd(R_tr,:)), 'Color', [1 0 0]);
b=plot(t_mids, mean(unch_ppd(R_tr,:)), 'Color', [0 0 1]);
legend([a,b],{'ch','unch'},'box','off','Location','northeast')

for i = 1:2
subplot(1,2,i)

    xlabel('time from pics (ms)')
ylabel('%\Delta post prob')

xlim([-500 1000])

title(['lever=',tt{i},' to chamber'])
end

set(gcf,'Position',[100 100 600 250])

% print([dataprefix,'_splitLR.png'],'-dpng')
end

