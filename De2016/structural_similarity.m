%A and B are cell arrays of atom neighbourhoods of molecules A and B.
function ss = structural_similarity(A,B)
gamma = 9;
zeta = 1;

n = length(A);
LSA = zeros(n,1);
for i = 1:n
    LSA(i) = local_similarity(A{i},A{i});
end

m = length(B);
LSB = zeros(m,1);
for j = 1:m
    LSB(j) = local_similarity(B{j},B{j});
end

C = zeros(n,m);
for i = 1:n
    for j = 1:m
        val = local_similarity(A{i},B{j})/sqrt(LSA(i)*LSB(j));
        C(i,j) = val^zeta;
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