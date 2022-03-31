function [unit_names] = add_unit_label(unit_names, chname,u)
% create new unit name, add to list

unit_letters = {'a','b','c','d','e','f'};

label = [chname,unit_letters{u}];

while sum(strcmp(label, unit_names))>0
    label = [label,unit_letters{u}];
end

unit_names(end+1,1) = {label};

end

