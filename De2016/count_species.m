%count max no of each atomic type in the dataset
n_mol = size(molecules,1)
max_types = zeros(5,1); %'H','C','N','O','S'
type2index = containers.Map({'H','C','N','O','S'},[1:5]);
for i = 1:n_mol
    types = char(molecules(i).atoms.types);
    count_types = zeros(5,1);
    for j = 1:size(types,1)
        count_types(type2index(types(j))) = count_types(type2index(types(j))) + 1;
    end
    max_types = max(max_types,count_types);
end
max_types
