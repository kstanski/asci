nt = 20;
train = repmat(Molecule(0), nt, 1);
ne = 10;
evaluate = repmat(Molecule(0), ne, 1);
[train,evaluate] = stratify(training_set_proper(50:(nt+ne)+49),train,evaluate);


lambda = 10^-15;
zeta = 1;
verbose = true;

tic
[~,RMSE,MAE,R2] = krr_de_radial(train,evaluate,lambda,zeta,verbose)
toc
