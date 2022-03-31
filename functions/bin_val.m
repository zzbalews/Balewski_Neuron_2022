function values_binned = bin_val(values,N)
% bin vector 'values' to 'N' bins (expect 16 unique items, 4 bins of 4 images)

possible = unique(values);
possible = possible(~isnan(possible));
nitems = length(possible)/N;

edges = [min(possible)-1,nan(1,N-1),max(possible)+1];
for i = 1:N-1
    edges(i+1) = mean(possible(nitems*i+[0,1]));
end

values_binned = discretize(values,edges);

end

