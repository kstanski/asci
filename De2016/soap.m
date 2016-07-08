%local similarity measure
function sum = soap(X1,X2)
l_max = 12;
sum = 0;
for l = 0:l_max
    sum = sum + I(l,X1,X2);
end
end

function s = I(l,X1,X2)
alpha = 0.4;

SH1 = spherical_harmonics(l,X1);
SH2 = spherical_harmonics(l,X2);

sx1 = size(X1,1);
sx2 = size(X2,1);
SB = zeros(sx1,sx2);
for idx1 = 1:sx1
    r1 = X1(idx1,:);
    for idx2 = 1:sx2
        r2 = X2(idx2,:);
        val = 4*pi*exp(-alpha*(norm(r1)^2+norm(r2)^2)/2);
        val = val * spherical_bessel(l,alpha*dot(r1,r2));
        SB(idx1,idx2) = val;
    end
end

s = 0;
for m1 = 1:2*l+1
    for m2 = 1:2*l+1
        sum = 0;
        for idx1 = 1:size(X1,1)
            for idx2 = 1:size(X2,1)
                v = SB(idx1,idx2);
                v = v*SH1(idx1,m1);
                v = v*conj(SH2(idx2,m2));
                sum = sum + v;
            end
        end
        s = s + conj(sum)*sum;
    end
end
end

function [theta,phi] = direction(r)
[theta,elev,~] = cart2sph(r(1),r(2),r(3));
phi = 2*pi - elev;
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

function SH = spherical_harmonics(l,X)
size_x = size(X,1);
SH = zeros(size_x,2*l+1);
for idx = 1:size_x
    [theta,phi] = direction(X(idx,:));
    lp = legendre(l,cos(theta));
    for m = -l:l
        if m > 0
            v = lp(m+1);
        else
            v = (-1)^m*factorial(l-m)/factorial(l+m)*lp(-m+1);
        end
        v = v * exp(sqrt(-1)*m*phi);
        SH(idx,m+l+1) = v * sqrt((2*l+1)/(4*pi)*factorial(l-m)/factorial(l+m));
    end
end
end
