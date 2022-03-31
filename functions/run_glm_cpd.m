function [ B_hat, B_se, coeff_p, model_p, CPD, SSR, SSR_const] = run_glm_cpd(X, Y)
%run_glm_cpd(X,Y)
% Input:
%   X = [1 x1 x2 ..] where x1, x2, ... are column explan vars
%   Y = column outcome var
% Output:
%   beta coeff estimates, etc.
% Y = X*B_hat+E, minimize E

N = size(X,1); % number trials

rk = rank(X); % rank
df = N - rk; % degrees of freedom

% find beta est
B_hat = pinv(X) * Y; % coeff est

% error
E = Y - X*B_hat; % error
SSR = diag(E'*E)'; % regression sum of squares

% find SE for beta est
s2 = sum(E.^2)/df; % var of population
B_se = diag(sqrt(inv(X'*X))).*sqrt(s2); % SE for coefficients

% find t-state, p-value for beta est
c = eye(size(B_hat,1)); % contrasts, if each coeff sig
temp = diag(c * pinv(X' * X) * c');
coeff_t = c*B_hat./sqrt(temp*s2);
coeff_p = 2*(1-tcdf(abs(coeff_t), df));

% compare to constant model
E_const = Y - mean(Y,1); % error
SSR_const = diag(E_const'*E_const)'; % regression sum of squares

nu_1 = rk-1;
nu_2 = N-rk;

F_stat = ((SSR_const - SSR)/nu_1./(SSR/nu_2));
model_p = 1-fcdf(F_stat,nu_1, nu_2);

% get CPD for each regressor-- compare to model w/o that var (skip
% procedure for constant, since did this above already)
CPD = nan(size(coeff_p));
for i = 2:size(X,2)
    X_subset = X;
    % use all regressors, except one being tested
    X_subset(:,i) = [];
    % get new beta est
    B_hat_subset = pinv(X_subset) * Y;
    % get new error
    E_subset = Y - X_subset*B_hat_subset;
    SSR_subset = diag(E_subset'*E_subset)';
    CPD(i,:) = (SSR_subset-SSR)./SSR_subset;
end

end

