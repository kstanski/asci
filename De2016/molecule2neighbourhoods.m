
function N = molecule2neighbourhoods(M)
coords = M.atoms.ff_coord;
n = size(coords,1);
dist = zeros(n);
threshold = 3;
for i = 1:n
    c_i = coords(i,:);
    for j = i:n
        c_j = coords(j,:);
        if norm(c_i-c_j) <= threshold
            dist(i,j) = 1;
            dist(j,i) = 1;
        end
    end
end

N = cell(n,1);
for i = 1:n
    nhood_i = zeros(sum(dist(i,:)),3);
    nhood_last = 0;
    for j = 1:n
        if dist(i,j) == 1
            nhood_last = nhood_last + 1;
            nhood_i(nhood_last,:) = coords(j,:) - coords(i,:);
        end
    end
    N{i} = nhood_i;
end
end
