function [output] = make_disc(x,edges)

mids = mean([edges(1:end-1);edges(2:end)]);

output = mids(discretize(x,edges));

end

