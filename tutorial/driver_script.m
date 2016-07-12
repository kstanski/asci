
kernel = @laplacian_kernel;
descriptor = @coulomb_matrix;
verbose = true;

lambda = 10^(-12);
sigma = 10^12;

nt = 10;
train = repmat(Molecule(0), nt, 1);
ne = 3;
evaluate = repmat(Molecule(0), ne, 1);
[train,evaluate] = stratify(hold_out_set(1:(nt+ne)),train,evaluate);
%[X,Y,X_p,Y_p] = apply_descriptor(train,evaluate,descriptor);

tic
[X,Y,X_p,Y_p] = apply_descriptor(training_set_proper,hold_out_set,descriptor);
[~,RMSE,MAE,R2] = krr(lambda,sigma,X,Y,X_p,Y_p,kernel,verbose)
toc
