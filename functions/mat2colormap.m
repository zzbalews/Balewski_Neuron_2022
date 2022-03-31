function X_smooth = mat2colormap(X)
X_smooth = [];
for i = 2:size(X,1)
    temp = [linspace(X(i-1,1),X(i,1),100)',...
        linspace(X(i-1,2),X(i,2),100)',...
        linspace(X(i-1,3),X(i,3),100)'];
    X_smooth = cat(1,X_smooth,temp);
end

end

