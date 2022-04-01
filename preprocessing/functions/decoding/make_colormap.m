function clrmp = make_colormap(X)

nsteps = size(X,1);

clrmp = [];
for s = 2:nsteps
    step = [];
    for i = 1:3
        temp = linspace(X(s-1,i),X(s,i),20)';
        step = cat(2,step,temp);
    end
    clrmp = cat(1,clrmp,step);
end
end
