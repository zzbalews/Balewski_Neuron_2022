function show_unit_proportions(unitcounts)

%%
subject_names = fieldnames(unitcounts);
region_names = fieldnames(unitcounts.(subject_names{1}));

nsubj = length(subject_names);
nreg = length(region_names);

varnames = {'val','dir','trtype'};
varnames_fancy = {'Max valuebin','Lever direction','Trial type'};
nvars = length(varnames);

clrs = [0 0 1; 1 1 0; 1 0 0]
figure;

for S = 1:nsubj
    
    subject = subject_names{S};
    disp(subject)
    subplot(1,2,S); hold on
    
    % plot prop sig units
    N = nan(nvars,nreg);
    prop = nan(size(N));
    Nsig = nan(1,nreg);
    
    for R = 1:nreg
        region = region_names{R};
        
        for v = 1:length(varnames)
            N(v,R) = unitcounts.(subject).(region).(varnames{v});
        end
        
        Nsig(R) = unitcounts.(subject).(region).Ntotal;
        
        prop(:,R) = N(:,R)/Nsig(R);
        
    end
    
    for v = 1:nvars
        b = bar(v,prop(v,:),'EdgeAlpha',0);
        b(1).FaceColor = clrs(v,:);
        b(2).FaceColor = clrs(v,:);
    end
    
    % add chi2 tests
    for v = 1:nvars
        observed = [ 1; -1] * N(v,:) + [0 0; Nsig];
        [chi2,df,chi2p] = do_chi2(observed);
        
        disp(['   ',varnames{v},': X2=',num2str(chi2),', df=',num2str(df),', p=',num2str(chi2p)]);
        if chi2p < 0.01/3
            
        text(v,0.55,'*','HorizontalAlignment','center')
        else
            text(v,0.55,'n.s.','HorizontalAlignment','center')
        end
    end
    
    
    % format
    
    set(gca,'XTick',1:3,'XTickLabel',varnames,'YTick',0:.25:.5)
    
    ylabel('proportion significant units')
    
    
end
        

end

