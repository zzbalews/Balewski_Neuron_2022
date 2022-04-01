function show_deltaRT_errortrials(bhvdata, tr_sampled, clrs)
% compare RT between error trials and correct trials, matched for offer
% value
%%
subject_names = fieldnames(tr_sampled);
trtypes = fieldnames(tr_sampled.(subject_names{1}));

figure;
Y = [0 1200; 0 750];

for s = 1:length(subject_names)
    
   subplot(1,2,s); hold on
   
   ll = [];
   
   rt1 = median(bhvdata.rt(tr_sampled.(subject_names{s}).(trtypes{1})));
   rt2 = median(bhvdata.rt(tr_sampled.(subject_names{s}).(trtypes{2})));
   
   rt_delta = rt1-rt2;
   disp([subject_names{s},': ',num2str(median(rt_delta)),'ms'])
%    for tr = 1:length(trtypes)
%        rt = median(bhvdata.rt(tr_sampled.(subject_names{s}).(trtypes{tr})));
       
       ci999 = quantile(rt_delta, [0.0005 0.5 0.9995]);
       ci999
       b=bar(1, ci999(2), 'FaceColor', clrs(1,:), 'EdgeAlpha', 0, 'FaceAlpha', .7);
       ll(1) = b;
       e = errorbar(1, ci999(2), ...
           diff(ci999(1:2)), ...
           diff(ci999(2:3)), 'k.');
%        if tr==length(trtypes)
           ll(end+1) = e;
%        end
           
%    end
%    ylim(Y(s,:))
%    
%    set(gca,'XTick',[1,2],'XTickLabel',{'correct','error'}, 'YTick', linspace(Y(s,1), Y(s,2), 4))
%       set(gca,'XTick',[],'YTick', linspace(Y(s,1), Y(s,2), 4))
% set(gca,'XTick',[],'YTick', linspace(Y(s,1), Y(s,2), 4))
   xlabel('free trials')
   ylabel('RT (ms)')
   title(subject_names{s})
   
   
%    legend(ll,{'correct','error','99% CI'},'box','off','Location','northwest')
   
end

end

