%A and B are arrays of atom neighbourhoods spectra of molecules A and B.
function ss = radial_structural_similarity(A,B,LSA,LSB)
gamma = 0.5;    %regularisation parameter

n = size(A,1);    %should be same as A
m = size(B,1);    %same as B

C = zeros(n,m);
for i = 1:n
    for j = 1:m
        val = radial_local_similarity(A(i,:),B(j,:))/sqrt(LSA(i)*LSB(j));
        C(i,j) = val;
    end
end
ss = trace(sinkhorn(C,gamma)'*C);
end


function P = sinkhorn(C,gamma)
K = exp((C-1)/gamma);
[n,m] = size(C);
en = repmat(1/n,n,1);
em = repmat(1/m,m,1);
v = em;
for r = 1:20
    u = en./(K*v);
    v = em./(K'*u);
end

P = zeros(n,m);
for i = 1:n
    for j = 1:m
        P(i,j) = u(i)*v(j)*K(i,j);
    end
end
end
