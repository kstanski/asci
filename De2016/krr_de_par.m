
function [f,RMSE,MAE,R2] = krr_de_par(training_set_proper,hold_out_set,lambda,verbose)
if verbose; disp('training'); end
if verbose; disp('computing neighbourhoods'); end
s = size(training_set_proper,1);
Ns = cell(s,1);
energy = zeros(s,1);
for idx = 1:s
    t = training_set_proper(s);
    Ns{idx} = molecule2neighbourhoods(t);
    energy(idx) = t.energy;
end

%Training
if verbose; disp('computing self similarity'); end
self_similarity = zeros(s,1);  %for efficient normalisation
parfor idx = 1:s
    self_similarity(idx) = structural_similarity(Ns{idx},Ns{idx});
end

if verbose; disp('computing structural similarity'); end
K = zeros(s);
for i = 1:s
    ns_elem = Ns{i};
    selfs_elem = self_similarity(i);
    parfor j = i+1:s
        k_ij = structural_similarity(ns_elem, Ns{j});
        k_ij = k_ij/sqrt(selfs_elem*self_similarity(j));
        K(i,j) = k_ij;
        %K(j,i) = k_ij;
    end
    K(i,i) = 1;
    K(i+1:s,i) = K(i,i+1:s)';
end

for i = 1:s
    ns_elem = Ns{i};
    selfs_elem = self_similarity(i);
    parfor j = 1:i
        k_ij = structural_similarity(ns_elem, Ns{j});
        k_ij = k_ij/sqrt(selfs_elem*self_similarity(j));
        K(i,j) = k_ij;
        %K(j,i) = k_ij;
    end
end

if verbose; disp('solving equation for alphas'); end
K = K + lambda*eye(s);
alpha = K\energy;


%Prediction
if verbose; disp('making prediction'); end
if verbose; disp('computing neighbourhoods'); end
s_p = size(hold_out_set,1);
Ns_p = cell(s_p,1);
energy_p = zeros(s_p,1);
for idx = 1:s_p
    h = hold_out_set(s_p);
    Ns_p{idx} = molecule2neighbourhoods(h);
    energy_p(idx) = h.energy;
end

if verbose; disp('computing self similarity'); end
self_similarity_p = zeros(s_p,1);  %for efficient normalisation
parfor idx = 1:s_p
    self_similarity_p(idx) = structural_similarity(Ns{idx},Ns{idx});
end

if verbose; disp('computing structural similarity'); end
L = zeros(s,s_p);
for i = 1:s
    ns_elem = Ns{i};
    selfs_elem = self_similarity(i);
    parfor j = 1:s_p
        l_ij = structural_similarity(ns_elem, Ns_p{j});
        l_ij = l_ij/sqrt(selfs_elem*self_similarity_p(j));
        L(i,j) = l_ij;
    end
end

f = L.'*alpha;

RMSE = rms(energy_p-f);
MAE = mean(abs(energy_p-f));
R2 = corr(energy_p,f)^2;

if verbose; plot_at_energy(energy_p,f); end

end
