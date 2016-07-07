%local similarity measure
function k = soap(X1,X2,zeta)
k = unnormalised_soap(X1,X2)/sqrt(unnormalised_soap(X1,X1)*unnormalised_soap(X2,X2));
k = k^zeta;
end

function sum = unnormalised_soap(X1,X2)
l_max = 2;
sum = 0;
for l = 0:l_max
    for m1 = 0:l
        for m2 = 0:l
            val = I(l,m1,m2,X1,X2);
            sum = sum + conj(val)*val;
        end
    end
end
end

function sum = I(l,m1,m2,X1,X2)
alpha = 0.4;
sum = 0;
for idx1 = 1:size(X1,1)
    r1 = X1(idx1,:);
    for idx2 = 1:size(X2,1)
        r2 = X2(idx2,:);
        sum = sum + I_t(l,m1,m2,alpha,r1,r2);
    end
end
end

function v = I_t(l,m1,m2,alpha,r1,r2)
v = 4*pi*exp(-alpha*(norm(r1)^2+norm(r2)^2)/2);
v = v*spherical_bessel(l,alpha*dot(r1,r2));
v = v*spherical_harmonic(l,m1,direction(r1));
v = v*conj(spherical_harmonic(l,m2,direction(r2)));
end

function r_unit = direction(r)
[theta,elev,~] = cart2sph(r(1),r(2),r(3));
phi = 2*pi - elev;
r_unit = [theta,phi];
end

%modified spherical bessel function of the first kind
function sb = spherical_bessel(l,x)
if abs(x) <= realmin
    if l == 0
        sb = 1;
    else
        sb = 0;
    end
else
    i = sqrt(-1);
    sb = sqrt(pi/(2*x)) * besselj(l+1/2,i*x) * i^-l;
end
end

function y = spherical_harmonic(l,m,r_unit)
theta = r_unit(1);
phi = r_unit(2);
y = legendre(l,cos(theta));
y = y(m+1);
y = y * exp(sqrt(-1)*m*phi);
y = y * sqrt((2*l+1)/(4*pi)*factorial(l-m)/factorial(l+m));
end
