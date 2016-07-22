%local similarity measure extended by atomic types
%S1 and S2 are atomic neighbourhoods spectra.
function sum = radial_local_similarity(S1,S2)
zeta = 1;
ntypes = length(S1);    %should be same for N2
K = eye(ntypes);
%K = diag([1,6,2,3,5]);
%K = ones(ntypes);    %alchemical similarity metric
%K(2:ntypes,1) = 0;    %'H','C','N','O','S'
%K(1,2:ntypes) = 0;
sum = 0;
for i = 1:ntypes
    for j = 1:ntypes
        if K(i,j) ~= 0
            sum = sum + K(i,j)*real(dot(S1{i},S2{j}))^zeta;
        end
    end
end
end
