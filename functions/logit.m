function output = logit(X)

X(X==1) = .99999;
X(X==0) = .00001;
output = log(X./(1-X));



end

