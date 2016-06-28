
kernel = @laplacian_kernel;
descriptor = @bag_of_bonds;

lambda = 10^(-6.5);
sigma = 724;

[RMSE,MAE,R2] = krr(lambda,sigma,training_set_proper,hold_out_set,kernel,descriptor)
