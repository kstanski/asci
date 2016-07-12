function g = radial_basis_function(r,r_cut,n)
n_max = 4;

S = zeros(n_max);
for i = 1:n_max
    for j = 1:n_max
        S(i,j) = sqrt((5+2*i)*(5+2*j))/(5+i+j);
    end
end

W = S^-0.5;

g = 0;
for alpha = 1:n_max
    g = g + W(n,alpha)*poly(r,r_cut,alpha);
end
end

function phi = poly(r,r_cut,alpha)
n_alpha = sqrt(r_cut^(2*alpha+5)/2*alpha+5);
phi = (r_cut-r)^(alpha+2)/n_alpha;
end
