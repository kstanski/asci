nt = 22;
train = repmat(Molecule(0), nt, 1);
ne = 10;
evaluate = repmat(Molecule(0), ne, 1);
[train,evaluate] = stratify(hold_out_set(1:(nt+ne)),train,evaluate);


lambda = 10^(-12);
zeta = 10^0;
verbose = true;

tic
[~,RMSE,MAE,R2] = krr_de(train,evaluate,lambda,zeta,verbose)
toc
