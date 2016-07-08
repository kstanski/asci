%local similarity measure extended by atomic types
%X1 and X2 are atomic neighbourhoods.
function sum = local_similarity(N1,N2)
ntypes = size(N1,1);    %should be same for N2
K = eye(ntypes);
%K = ones(ntypes);    %alchemical similarity metric
%K(2:ntypes,1) = 0;    %'H','C','N','O','S'
%K(1,2:ntypes) = 0;
sum = 0;
for i = 1:ntypes
    for j = 1:ntypes
        if K(i,j) ~= 0
            sum = sum + K(i,j)*soap(N1{i},N2{j});
        end
    end
end
end
