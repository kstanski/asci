
function [f,RMSE,MAE,R2] = krr_de(training_set_proper,hold_out_set,lambda,zeta,verbose)
if verbose; disp('training'); end
if verbose; disp('computing neighbourhoods'); end
s = size(training_set_proper,1);
Ns = cell(s,1);
energy = zeros(s,1);
for idx = 1:s
    t = training_set_proper(idx);
    Ns{idx} = molecule2neighbourhoods(t);
    energy(idx) = t.energy;
end

%Training
if verbose; disp('computing self similarity'); end
self_similarity = zeros(s,1);  %for efficient normalisation
for idx = 1:s
    nhood = Ns{idx};
    self_similarity(idx) = structural_similarity(nhood,nhood);
end

if verbose; disp('computing structural similarity'); end
K = zeros(s);
for i = 1:s
    disp(i);
    K(i,i) = 1;
    for j = i+1:s
        k_ij = structural_similarity(Ns{i}, Ns{j});
        k_ij = k_ij/sqrt(self_similarity(i)*self_similarity(j));
        K(i,j) = k_ij;
        K(j,i) = k_ij;
    end
end

if verbose; disp('solving equation for alphas'); end
K = K.^zeta + lambda*eye(s)
alpha = K\energy;


%Prediction
if verbose; disp('making prediction'); end
if verbose; disp('computing neighbourhoods'); end
s_p = size(hold_out_set,1);
Ns_p = cell(s_p,1);
energy_p = zeros(s_p,1);
for idx = 1:s_p
    h = hold_out_set(idx);
    Ns_p{idx} = molecule2neighbourhoods(h);
    energy_p(idx) = h.energy;
end

if verbose; disp('computing self similarity'); end
self_similarity_p = zeros(s_p,1);  %for efficient normalisation
for idx = 1:s_p
    nhood_p = Ns_p{idx};
    self_similarity_p(idx) = structural_similarity(nhood_p,nhood_p);
end

if verbose; disp('computing structural similarity'); end
L = zeros(s,s_p);
for i = 1:s
    disp(i);
    for j = 1:s_p
        l_ij = structural_similarity(Ns{i}, Ns_p{j});
        l_ij = l_ij/sqrt(self_similarity(i)*self_similarity_p(j));
        L(i,j) = l_ij;
    end
end

L = L.^zeta;
f = L'*alpha;

RMSE = rms(energy_p-f);
MAE = mean(abs(energy_p-f));
R2 = corr(energy_p,f)^2;

if verbose; plot_at_energy(energy_p,f); end

end
