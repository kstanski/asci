
function k = laplacian_kernel(x_i,x_j,sigma)
    k = exp(-norm(x_i-x_j,1)/sigma);
end
