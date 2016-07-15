nt = 15;
train = repmat(Molecule(0), nt, 1);
ne = 5;
evaluate = repmat(Molecule(0), ne, 1);
[train,evaluate] = stratify(hold_out_set(1:(nt+ne)),train,evaluate);


lambda = 10^-15;
zeta = 1;
verbose = true;

tic
[~,RMSE,MAE,R2] = krr_de_radial(train,evaluate,lambda,zeta,verbose)
toc
