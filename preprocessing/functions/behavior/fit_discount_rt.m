function [w_names,w_fit,BIC,R2] = ...
    fit_discount_rt(functype,lever,jpg,amnt,prob,rt)

%% reduce to unique trial rows
% by amnt/prob offer pair + response time
TR = struct();

temp = unique([jpg,rt],'rows');
TR.jpg = temp(:,1:2);
TR.rt = temp(:,3);

for i = 1:size(TR.jpg)

    % get matching trial rows
    idx = sum(jpg==TR.jpg(i,:),2)==2 & ...
        rt==TR.rt(i);
    
    % amnt/prob for trials
    TR.amnt(i,:) = amnt(find(idx,1),:);
    TR.prob(i,:) = prob(find(idx,1),:);
    
    % get # left/right lever
    idx_leftlever = idx & lever==-1;
    idx_rightlever = idx & lever==+1;
    TR.N(i,:) = [sum(idx_leftlever),sum(idx_rightlever)];
    
end

%% do fit

ntr = length(lever);

if strcmp(functype,'lin') % linear
w_names = {'A - w1(1-P)'};
elseif strcmp(functype,'hyp') % hyperbolic
w_names = {'A / 1 + w1(1-P)'};
elseif strcmp(functype,'exp') % exponential 
w_names = {'A * exp(w1(1-P))'};
elseif strcmp(functype,'expval') % expected value
w_names = {};
end

w_names(end+1:end+3) = {'softmax w2*deltaval','softmax w3*rt','softmax side bias'};

loglikefunc = @(w) ...
    -sum(log(softmax_discount_rt(functype,{TR.amnt,TR.prob},TR.rt,w) .^ TR.N(:,1))) + ...
    -sum(log((1-softmax_discount_rt(functype,{TR.amnt,TR.prob},TR.rt,w)) .^ TR.N(:,2)));


[w_fit,ll_fit] = best_fminsearch(loglikefunc,w_names);

BIC = length(w_fit)*log(ntr) +2*ll_fit;

%% get fit R2
valuefunction = eval(['@',functype,'discount']);
if strcmp(functype,'expval')
subjval = valuefunction(amnt,prob);
else
subjval = valuefunction(amnt,prob,w_fit(1));
end

V_delta = subjval(:,1) - subjval(:,2);

goleft = lever==-1; % 1=left,0=right
goleft_avg = mean(goleft); % prop trials left lever overall
goleft_fit = softmax_discount_rt(functype,{V_delta},rt,w_fit);

R2 = 1 - sum((goleft-goleft_fit).^2) / sum((goleft-goleft_avg).^2);

end

