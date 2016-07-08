

%tsp = training_set_proper(1);
%tsp.atoms.ff_coord = tsp.atoms.ff_coord + 1;

N1 = molecule2neighbourhoods(training_set_proper(5));
%N1 = N1{1};
N2 = molecule2neighbourhoods(training_set_proper(5));
%N2 = N2{1};

tic
lso = zeros(size(N1,1),size(N2,1));
for i = 1:size(N1,1)
    n1 = N1{i};
    for j = 1:size(N2,1)
        n2 = N2{j};
        lso(i,j) = local_similarity_old(n1,n2);
    end
end
toc

tic
ls = zeros(size(N1,1),size(N2,1));
for i = 1:size(N1,1)
    n1 = N1{i};
    for j = 1:size(N2,1)
        n2 = N2{j};
        ls(i,j) = local_similarity(n1,n2);
    end
end
toc

for i = 1:size(N1,1)
    for j = 1:size(N2,1)
        diff = lso(i,j) - ls(i,j);
        if abs(diff) > 10^-15 || isnan(diff)
            disp(diff)
        end
    end
end

