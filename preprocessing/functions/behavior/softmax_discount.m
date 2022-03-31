function output = softmax_discount(type,input,w)

valuefunction = eval(['@',type,'discount']);

if length(input)==2 % input is amnt and prob
    
    amnt = input{1};
    prob = input{2};
    
    if strcmp(type,'expval')
    V_left = valuefunction(amnt(:,1),prob(:,1));
    V_right = valuefunction(amnt(:,2),prob(:,2));
    else
    V_left = valuefunction(amnt(:,1),prob(:,1),w(1));
    V_right = valuefunction(amnt(:,2),prob(:,2),w(1));
    end
    
    V_delta = V_left-V_right;
    
elseif length(input)==1 % input is subj val difference
    
    V_delta = input{1};
    
end

output = 1./(1 + exp( w(end-1)*V_delta + w(end) ));

end

