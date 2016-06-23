
function k = gaussian_kernel(x_i,x_j,sigma)
    k = exp(-(norm(x_i-x_j)^2)/(2*sigma^2));
end
