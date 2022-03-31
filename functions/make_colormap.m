function clrmap = make_colormap(X)

nsteps = size(X,1);

clrmap = [];
for s = 2:nsteps
    step = [];
    for i = 1:3
        temp = linspace(X(s-1,i),X(s,i),20);
        step = cat(2,step,temp');
    end
    clrmap = cat(1,clrmap,step);
end
end
