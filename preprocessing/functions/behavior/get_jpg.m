function jpg = get_jpg(amnt,prob)

amnt(amnt==0) = NaN;
prob(prob==0) = NaN;

amnt_disc = [];
prob_disc = [];

for i = 1:2
    
[~,~,temp] = unique(amnt(:,i));
temp(temp>4) = 5;
amnt_disc(:,i) = temp;

[~,~,temp] = unique(prob(:,i));
temp(temp>4) = 5;
prob_disc(:,i) = temp;

end

jpg_ID = [...
    1, 2, 3, 4, 0;
    5, 6, 7, 8, 0;
    9,10,11,12, 0;
    13,14,15,16,0;
    0, 0, 0, 0, 17];

jpg = [diag(jpg_ID(amnt_disc(:,1),prob_disc(:,1))), ...
    diag(jpg_ID(amnt_disc(:,2),prob_disc(:,2)))];
end

