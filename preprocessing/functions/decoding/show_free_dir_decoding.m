function show_free_dir_decoding(decode_output,minidata,decoder)

%% load dir decoder outputs

load([decode_output,'.mat'])

nmids = length(t_mids);
%% sort free trials by RT
if strcmp(decoder,'lever')
truelabel = minidata.lever;
elseif strcmp(decoder,'sacc1')
    truelabel = minidata.saccloc(:,1);
    elseif strcmp(decoder,'sacc2')
    truelabel = minidata.saccloc(:,2);
end

free = find(minidata.trialtype==2 & minidata.lever~=0 & isnan(truelabel));
[rt_sort,idx_sort] = sort(minidata.rt(free));

free_sort = free(idx_sort);

ntr = length(free_sort);

%% rearrange by chosen/unchosen

postprob_chosen = alltrials_postprob_avg(:,:,:,1);
postprob_chosen(truelabel==1,:) = alltrials_postprob_avg(truelabel==1,:,:,2);

state_chosen = alltrials_states(:,:,:,1,5);
state_chosen(truelabel==1,:) = alltrials_states(truelabel==1,:,:,2,5);

state_unchosen = alltrials_states(:,:,:,2,5);
state_unchosen(truelabel==1,:) = alltrials_states(truelabel==1,:,:,1,5);

%% show!


rb = redbluecmap;
rb = rb+.2;
rb(rb>1)= 1;

for i = 1:3
    
    figure;
    if i==1
        imagesc(postprob_chosen(free_sort,:))
        colormap(rb);
        title('post. prob. for chosen direction')
        
    elseif i==2
        imagesc(state_chosen(free_sort,:))
        heat = make_colormap([1 1 1;1 0 0 ]);
        colormap(heat);
        title('chosen dir states')
        
    elseif i==3
        imagesc(state_unchosen(free_sort,:))
        heat = make_colormap([1 1 1;0 0 1]);
        colormap(heat);
        title('alt dir states')
    end
    hold on
    plot(399+rt_sort/5,1:ntr,'k','LineWidth',2)
    
    plot([399 399],0.5+[0,ntr],'k')
    set(gca,'XTick',49:50:797,'XTickLabel',1+t_mids(49:50:797))
    
    xlim([298,651])
    
    xlabel('time from pics (ms)')
    ylabel('trials (order by RT)')
    set(gcf,'Position',[100 100 550 1100])
end
end


