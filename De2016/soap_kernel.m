
function k = soap_kernel(X1,X2)

end

%X environment (positions of atoms)
%r position vector
function sum = local_density(X,r,sigma)
    sum = 0;
    for x = X'
        sum = sum + gaussian(x,r,sigma);
    end
end


%x center of atom
%r position vector
function g = gaussian(x,r,sigma)
    g = exp(-norm(x-r)^2/(2*sigma^2));
end
