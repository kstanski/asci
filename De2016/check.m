
tsp = training_set_proper(1);
n_at = size(tsp.atoms.ff_coord,1);
%permutation
%tsp.atoms.ff_coord = vertcat(tsp.atoms.ff_coord(1:2,:),tsp.atoms.ff_coord(3:n_at,:));
%translation
%tsp.atoms.ff_coord = tsp.atoms.ff_coord + 1.5;
%rotation
tsp.atoms.ff_coord = quatrotate([pi/0.9,pi/2.2,pi/3.01,pi/4.1],tsp.atoms.ff_coord);

%arbitrary alteration
%tsp.atoms.ff_coord(2,3) = tsp.atoms.ff_coord(2,3) + 0.1;

N1 = molecule2neighbourhoods(training_set_proper(1));
N2 = molecule2neighbourhoods(tsp);

%ss = structural_similarity(N1,N2);
%ss = ss/sqrt(structural_similarity(N1,N1)*structural_similarity(N2,N2));
%disp(sprintf('%.10f',ss))

lambda = 10^-10;
verbose = true;

%train = hold_out_set(1:2);
%train = vertcat(train,hold_out_set(4:5));
%t = molecule2neighbourhoods(train);

%evaluate = hold_out_set(3);
%evaluate = vertcat(evaluate,hold_out_set(5));
%e = molecule2neighbourhoods(evaluate);

%structural_similarity(t,e)/sqrt(structural_similarity(t,t)*structural_similarity(e,e))
train = repmat(Molecule(0), 50, 1);
evaluate = repmat(Molecule(0), 20, 1);
[train,evaluate] = stratify(hold_out_set(1:70),train,evaluate);

tic
[~,RMSE,MAE,R2] = krr_de(train,evaluate,lambda,verbose)
toc

%tic
%[f,~,~,~] = krr_de_par(hold_out_set(3:5),hold_out_set(6),lambda,verbose)
%toc