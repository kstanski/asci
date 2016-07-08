%N is a collection of neighbourhoods.
%Every neighbourhood is a cell array,
%in which each cell corresponds to an atomic types and holds a matrix of
%coordinates of neighbouring atoms of this particular type. 
function N = molecule2neighbourhoods(M)

max_types_counts = [16,7,3,3,1]';    %depends on the dataset used

ntypes = 5;
type2index = containers.Map({'H','C','N','O','S'},[1:ntypes]);
types = char(M.atoms.types);   %'H','C','N','O','S'
coords = M.atoms.ff_coord;
centers = sum(max_types_counts(1:ntypes));

n = size(coords,1);
neighbour_matrix = zeros(n);
no_each_type = zeros(n,ntypes);
cutoff = 3;
for i = 1:n
    c_i = coords(i,:);
    for j = i:n
        c_j = coords(j,:);
        if norm(c_i-c_j) <= cutoff
            neighbour_matrix(i,j) = 1;
            no_each_type(i,type2index(types(j))) = no_each_type(i,type2index(types(j))) + 1;
            if i ~= j
                neighbour_matrix(j,i) = 1;
                no_each_type(j,type2index(types(i))) = no_each_type(j,type2index(types(i))) + 1;
            end
        end
    end
end

N = cell(centers,1);
nlast = 0;
types_count = zeros(ntypes,1);
for i = 1:n
    t = type2index(types(i));
    types_count(t) = types_count(t) + 1;
    if t > 0    %include all species as centers
        nhood_i = cell(ntypes,1);
        for type = 1:ntypes
            nhood_i{type} = zeros(no_each_type(i,type),3);
        end
        nhood_last = zeros(ntypes,1);
        for j = 1:n
            if neighbour_matrix(i,j) == 1
                type_idx = type2index(types(j));
                nhood_last(type_idx) = nhood_last(type_idx) + 1;
                nhood_i{type_idx}(nhood_last(type_idx),:) = coords(j,:) - coords(i,:);
            end
        end
        nlast = nlast + 1;
        N{nlast} = nhood_i;
    end
end

top_up = max_types_counts - types_count;
for t = 1:ntypes
    for i = 1:top_up(t)
        nlast = nlast + 1;
        nhood_i = cell(ntypes,1);
        nhood_i{t} = [0,0,0];
        N{nlast} = nhood_i;
    end
end

end
