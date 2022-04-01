function unitcounts = show_unit_venn(glmdata, t_period, out_prefix)

region_names = unique(glmdata.region);
subject_names = unique(glmdata.subject);

nreg = length(region_names);
nsubj = length(subject_names);

useful_time = glmdata.t_mids>=t_period(1) & glmdata.t_mids<=t_period(2);
useful_vars = glmdata.glmvars(2:4);

textout = fopen([out_prefix,'.txt'],'w'); % save counts in text file, as backup

label = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; 0 0 0];

unitcounts = struct();

figure;
for R = 1:nreg
    for S = 1:nsubj
        
        % restrict units
        region = region_names{R};
        subject = subject_names{S};
                fprintf(textout,'\n\n%s %s\n',region,subject);

        idx = strcmp(glmdata.region,region_names{R}) & ...
            strcmp(glmdata.subject,subject_names{S});
        
        % look at 't_period' chunk (ignore constant & trialnum vars)
        subset = glmdata.sigunit(idx,useful_time,2:4);
        ifsig = permute(sum(subset,2)>0,[1 3 2]);

        % get unique counts
        label_count = nan(8,1);
        for L = 1:8
            label_count(L) = sum(sum(ifsig==label(L,:),2)==3);
            
            % update text file (backup)
            temp = cellfun(@(x) [x,', '],useful_vars(label(L,:)==1),'UniformOutput',false);
            fprintf(textout,'%d %s\n',label_count(L),[temp{:}]);
            
        end
          
        fprintf(textout,'%s %d%s%d %s\n','total sig: ',sum(sum(ifsig,2)>0),'/',sum(idx),...
            ['(',num2str(round(100*sum(sum(ifsig,2)>0)/sum(idx))),'%)']);
        
        % update unit_counts
        unitcounts.(subject).(region).val = sum(label_count(label(:,3)==1));
        unitcounts.(subject).(region).dir = sum(label_count(label(:,2)==1));
        unitcounts.(subject).(region).trtype = sum(label_count(label(:,1)==1));
        unitcounts.(subject).(region).Nsig = sum(label_count(1:7));
        unitcounts.(subject).(region).Ntotal = sum(label_count);
        
        % plot
        subplot(nreg,nsubj,S + (R-1)*nsubj); hold on
        [~,venninfo] = venn(label_count(1:7));
                
        % add text
        for L = 1:7
            text(venninfo.ZoneCentroid(L,1),...
                venninfo.ZoneCentroid(L,2),...
                num2str(label_count(L)), 'HorizontalAlignment','center');
        end
        
        % format
        daspect([1 1 1])
        title([region,' ', subject])
        
        % do chi2 test: are more lever&value than expected by chance?
        nunits_val = sum(label_count(label(:,3)==1));
        nunits_dir = sum(label_count(label(:,2)==1));
        nunits_valdir = sum(label_count(label(:,2)==1 & label(:,3)==1));
        nunits_neither = sum(label_count(label(:,2)~=1 & label(:,3)~=1));
        
        O = [nunits_valdir, nunits_val-nunits_valdir;
            nunits_dir-nunits_valdir, nunits_neither];
        E = repmat(sum(O),2,1) .* repmat(sum(O,2),1,2)/sum(sum(O));
        chi2 = sum(sum(((O-E).^2)./E));
        df = (2-1)*(2-1);
        chi2p = 1 - chi2cdf(chi2,df);
        
        disp([subject,' @ ',region,': # val & dir units chi2 test'])
        disp(['... chi2=',num2str(chi2),', df=1, p=',num2str(chi2p)])
        
    end
end

% save fig
print([out_prefix,'.svg'],'-dsvg')

% close backup text file
fclose(textout)


end

