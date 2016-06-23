function k = test_kernel(x_i,x_j,sigma)
    p = 1.326;
    k = exp(-(norm(x_i-x_j,p)^p)/(p*sigma^p));
end