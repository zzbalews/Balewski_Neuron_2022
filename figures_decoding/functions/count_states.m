function [t_idx,track_flips]=count_states(bhvdata,decdata,keep_tr, decoder)
%%

thresh = 0.01;

subject_names = unique(decdata.subject);
nsubj = length(subject_names);

if strcmp(decoder,'value')
    STATES_ch = decdata.ch_state;
    STATES_unch = decdata.unch_state;
    STATES_na = decdata.na_state;
    flipLIM = [0 10];
    optLIM = [0 4];
    optLABEL = {'ch','unch','na'};
    terms = 3;
elseif strcmp(decoder,'direction')
   
    [ntr,nmids] = size(decdata.postprob_ch_disc);
    STATES_ch = decdata.postprob_ch_disc==1;
    STATES_unch = decdata.postprob_ch_disc==-1;
    STATES_na = nan(ntr,nmids,2);
    flipLIM = [0 5];
    optLIM = [0 3];
    optLABEL = {'ch','unch'};
    terms = 2;
    
end

track_flips = nan(size(STATES_ch,1),2);


% useful vars
t_mids = decdata.t_mids;

figure;
T = 800;

t_idx = find(t_mids>=0 & t_mids<T);
for s = 1:nsubj
    
    % restrict trials
    subject = subject_names{s};
    disp([subject,' value decoding'])
    tr_idx = find(strcmp(bhvdata.subject,subject) & ...
        keep_tr);
    
    %% count states
    ntr = length(tr_idx);
    count_ch = nan(ntr,1);
    dur_ch = nan(ntr,1);
    
    count_unch = nan(ntr,1);
    dur_unch = nan(ntr,1);
    
    count_na = nan(ntr,1);
    dur_na = nan(ntr,1);
    
    temp = repelem([1,2],floor(ntr/2),ceil(ntr/2));
    which_na = temp(randperm(ntr));
    
    count_flips = nan(ntr,1);
    
    track_states = cell(ntr,1);
    for tr = 1:ntr
        states = [];
        
        % ch
        temp = diff([0, STATES_ch(tr_idx(tr),t_idx),0]);
        ch = find(temp==1);
        ch_off = find(temp==-1);
        
        count_ch(tr) = length(ch);
        if length(ch)>0
            states = cat(1,states,[repmat(1,length(ch),1),ch']);
            dur_ch(tr) = mean((ch_off-ch)*5+15);
        end
        
        % unch
        temp = diff([0, STATES_unch(tr_idx(tr),t_idx),0]);
        unch = find(temp==1);
        unch_off = find(temp==-1);
        
        count_unch(tr) = length(unch);
        if length(unch)>0
            states = cat(1,states,[repmat(-1,length(unch),1),unch']);
            dur_unch(tr) = mean((unch_off-unch)*5+15);
        end
        
        % na
        temp = diff([0, STATES_na(tr_idx(tr),t_idx,which_na(tr)),0]);
        na = find(temp==1);
        na_off = find(temp==-1);
        
        count_na(tr) = length(na);
        if length(na)>0
            %             states = cat(1,states,[repmat(-1,length(na),1),na']);
            dur_na(tr) = mean((na_off-na)*5+15);
        end
        
        
        if size(states,1)>0
            [~,idx] = sort(states(:,2));
            states = states(idx,:);
            track_states{tr} = states;
            
            check_flips = [0; diff(states(:,1))];
            
            count_flips(tr) = sum(check_flips~=0);
            if sum(check_flips~=0>0)
                temp = find(check_flips~=0,1);
                track_flips(tr_idx(tr),:) = [states(temp,2), check_flips(temp)];
           
            end
        end
    end
    disp(['% trials with at least 1 ch/unch state:'])
    disp(['   ',subject,':',num2str(round(100*mean(~isnan(count_flips)),1)),'%'])
    %% show # flips
    subplot(1,6,4+s);
    histogram(count_flips,'BinWidth',1,'FaceColor',[.7 0 .7],'EdgeAlpha',0,'Normalization','probability');
    disp(['   median # flips (ch & unch): ',num2str(nanmedian(count_flips))])
    xlabel('# flips between ch & unch')
    ylabel('norm. trial count')
    xlim(flipLIM)
%     ylim([0 .8])
ylim([0 .4])
%     set(gca,'box','off','YTick',0:.4:.8)
    set(gca,'box','off','YTick',0:.2:.4)
    
    
   %% plot # of states
   disp('# states/trial:')
   temp = count_ch+count_unch+count_na;
   disp(['   ',subject,': ',num2str(round(mean(temp),1)),' +- ',num2str(round(std(temp)/sqrt(length(temp)),2))])
    clrs_tr = [1 0 0; 0 0 1; .7 .7 .7];
    
    subplot(1,6,s); hold on
    temp = [count_ch,count_unch,count_na];
    temp = temp(:,1:terms);
    for i = 1:terms
        bar(i,mean(temp(:,i)),'FaceColor',clrs_tr(i,:),'EdgeAlpha',0);
    end
    errorbar(1:terms,mean(temp),std(temp)/sqrt(ntr),'k.')
    set(gca,'XTick',1:3,'XTickLabel',{'Ch','UnCh','NA'})
    
    % # trials with ch, unch, both
    temp = [count_ch,count_unch]>0;
    
    Nch = sum(sum(temp==[1 0],2)==2);
    Nunch = sum(sum(temp==[0 1],2)==2);
    Nboth = sum(sum(temp==[1 1],2)==2);
    Nnone = sum(sum(temp==[0 0],2)==2);
    
    O = [Nboth, Nch; Nunch, Nnone];
    E = repmat(sum(O),2,1) .* repmat(sum(O,2),1,2)/sum(sum(O));
    chi2 = sum(sum(((O-E).^2)./E));
    df = (2-1)*(2-1);
    chi2p = 1 - chi2cdf(chi2,df);
    
    disp('   more ch & unch than expected by chance?:')
    disp(['   X2=',num2str(chi2),', df=',num2str(df),', p=',num2str(chi2p)]);
    
    % anova
    data = table();
    data.count = [count_ch;count_unch;count_na];
    data.option = repelem({'ch';'unch';'other'},ntr);
    data.tr = repmat([1:ntr]',3,1);
    
    
    [P,ANOVATAB,STATS] = anova1(data.count, data.option,'off');
    disp('   # states:')
    disp(['   ANOVA: F(',num2str(ANOVATAB{2,3}),',',num2str(ANOVATAB{3,3}),')=',num2str(round(ANOVATAB{2,5},1)),...
        ',p<1e',num2str(ceil(log10(P)))])
    
    Y = 2.75;
    
         [H,P,CI,STATS] = ttest(count_ch,count_unch);
         if P<thresh
             plot([1.1 1.9],[Y Y],'k')
             text(1.5,Y+.25,['p<1e',num2str(ceil(log10(P)))],'HorizontalAlignment','left','Rotation',90)
         end
         
  if terms==3
        [H,P,CI,STATS] = ttest(count_unch,count_na);
        if P<thresh
             plot([2.1 2.9],[Y Y],'k')
             text(2.5,Y+.25,['p<1e',num2str(ceil(log10(P)))],'HorizontalAlignment','left','Rotation',90)
        end
  end
    % format
    set(gca,'XTick',1:terms,'XTickLabel',optLABEL,'YTick',1:5)
    xlim(optLIM)
%     ylim([0 3])
 ylim([0 4.6])
    ylabel(['# states in ',num2str(T),' ms pics window'])
    
    
    
    %% plot state durs
     subplot(1,6,2+s); hold on
    temp = [dur_ch,dur_unch,dur_na];
    temp = temp(:,1:terms);
    for i = 1:terms
        bar(i,nanmean(temp(:,i)),'FaceColor',clrs_tr(i,:),'EdgeAlpha',0);
    end
    errorbar(1:terms,nanmean(temp),nanstd(temp)./sqrt(sum(~isnan(temp))),'k.')
    plot([0 4],(15+4*5)*[1 1],'k:')
   
        % anova
    data = table();
    data.dur = [dur_ch;dur_unch;dur_na];
    data.option = repelem({'ch';'unch';'other'},ntr);
    data.tr = repmat([1:ntr]',3,1);
    
    [P,ANOVATAB,STATS] = anova1(data.dur, data.option,'off');
    disp('   state duration:')
    disp(['   ANOVA: F(',num2str(ANOVATAB{2,3}),',',num2str(ANOVATAB{3,3}),')=',num2str(round(ANOVATAB{2,5},1)),...
        ',p<1e',num2str(ceil(log10(P)))])
    
    Y = 95;
             [H,P,CI,STATS] = ttest(dur_ch,dur_unch);
         if P<thresh
             plot([1.1 1.9],[Y Y],'k')
             text(1.5,Y+5,['p<1e',num2str(ceil(log10(P)))],'HorizontalAlignment','left','Rotation',90)
         end
         
         if terms==3
  
        [H,P,CI,STATS] = ttest(dur_unch,dur_na);
        if P<thresh
             plot([2.1 2.9],[Y Y],'k')
             text(2.5,Y+5,['p<1e',num2str(ceil(log10(P)))],'HorizontalAlignment','left','Rotation',90)
        end
         end
    set(gca,'XTick',1:terms,'XTickLabel',optLABEL,'YTick',0:50:100)
    xlim(optLIM)
%     ylim([0 110])
     ylim([0 85])
    ylabel('state dur (ms)')
    
end
%%


end