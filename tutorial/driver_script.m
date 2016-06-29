
kernel = @laplacian_kernel;
descriptor = @bag_of_bonds;
verbose = true;

lambda = 10^(-13);
sigma = 10^11;

[X,Y,X_p,Y_p] = apply_descriptor(training_set_proper,hold_out_set,descriptor);
[~,RMSE,MAE,R2] = krr(lambda,sigma,X,Y,X_p,Y_p,kernel,verbose)
