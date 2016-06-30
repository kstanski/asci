
kernel = @laplacian_kernel;
descriptor = @bag_of_bonds;
verbose = true;

lambda = 10^(-12);
sigma = 10^12;

[X,Y,X_p,Y_p] = apply_descriptor(training_set_proper,hold_out_set,descriptor);
[~,RMSE,MAE,R2] = krr(lambda,sigma,X,Y,X_p,Y_p,kernel,verbose)
