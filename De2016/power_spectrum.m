%rotationally invariant power spectrum
function P = power_spectrum(X)

l_max = 10;      %spherical harmonics
n_max = 9;      %radial basis functions
r_cut = 3;

[Theta,Phi,R] = spherical_coords(X);

RBF = zeros(length(R),n_max);
for idx = 1:length(R)
    for n = 1:n_max
        RBF(idx,n) = radial_basis_function(R(idx),r_cut,n,n_max);
    end
end

P = zeros(l_max+1,n_max,n_max);
for l = 0:l_max
    SH = spherical_harmonics(l,Theta,Phi);
    for n1 = 1:n_max
        for n2 = 1:n1
            for m = 1:2*l+1
                c1 = 0;
                c2 = 0;
                for idx = 1:length(R)
                    c1 = c1 + RBF(idx,n1)*SH(idx,m);
                    c2 = c2 + RBF(idx,n2)*SH(idx,m);
                end
                P(l+1,n1,n2) = P(l+1,n1,n2) + conj(c1)*c2;
                if n1 ~= n2
                    P(l+1,n2,n1) = P(l+1,n2,n1) + conj(c2)*c1;
                end
            end
        end
    end
end

P = P(:);
if norm(P) ~= 0
    P = P/norm(P);
end
P = real(P);
end

function [Theta,Phi,R] = spherical_coords(X)
if size(X,2) == 3
    [Theta,Elev,R] = cart2sph(X(:,1),X(:,2),X(:,3));
else
    Theta = 0;
    Elev = 0;
    R = 0;
end
Phi = 2*pi - Elev;
end

function SH = spherical_harmonics(l,Theta,Phi)
size_x = size(Theta,1);    %should be same as Phi
SH = zeros(size_x,2*l+1);
for idx = 1:size_x
    theta = Theta(idx);
    phi = Phi(idx);
    lp = legendre(l,cos(theta));
    for m = 0:l
        v = lp(m+1);
        v = v * exp(sqrt(-1)*m*phi);
        v = v * sqrt((2*l+1)/(4*pi)*factorial(l-m)/factorial(l+m));
        SH(idx,l+1+m) = v;
        if m > 0
            SH(idx,l+1-m) = -1^-m * conj(v);
        end
    end
end
end

function g = radial_basis_function(r,r_cut,n,n_max)
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
