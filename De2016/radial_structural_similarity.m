%A and B are arrays of atom neighbourhoods spectra of molecules A and B.
function ss = radial_structural_similarity(A,B)
gamma = 0.5;    %regularisation parameter

n = size(A,1);
LSA = zeros(n,1);
for i = 1:n
    LSA(i) = radial_local_similarity(A(i,:),A(i,:));
end

m = size(B,1);
LSB = zeros(m,1);
for j = 1:m
    LSB(j) = radial_local_similarity(B(j,:),B(j,:));
end

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
