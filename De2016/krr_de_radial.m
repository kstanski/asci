
function [f,RMSE,MAE,R2] = krr_de_radial(training_set_proper,hold_out_set,lambda,zeta,verbose)
if verbose; disp('training'); end

if verbose; disp('computing spectra'); end
s = size(training_set_proper,1);
Ss = cell(s,1);    %spectras
energy = zeros(s,1);
parfor idx = 1:s
    %disp(idx);
    t = training_set_proper(idx);
    N = molecule2neighbourhoods(t);
    Ss{idx} = neighbourhoods2spectra(N);
    energy(idx) = t.energy;
end

if verbose; disp('computing local similarity'); end
LSs = cell(s,1);    %local similarities
parfor idx = 1:s
    spectra = Ss{idx};
    n = size(spectra,1);
    LS = zeros(n,1);
    for i = 1:n
        LS(i) = radial_local_similarity(spectra(i,:),spectra(i,:));
    end
    LSs{idx} = LS;
end

if verbose; disp('computing self structural similarity'); end
self_similarity = zeros(s,1);  %for efficient normalisation
parfor idx = 1:s
    spectra = Ss{idx};
    LS = LSs{idx};
    self_similarity(idx) = radial_structural_similarity(spectra,spectra,LS,LS);
end

if verbose; disp('computing cross structural similarity'); end

K = zeros(s);
parfor i = 1:s
    %K(i,i) = 1;
    for j = 1:s
        if j > i
            k_ij = radial_structural_similarity(Ss{i}, Ss{j}, LSs{i}, LSs{j});
            k_ij = k_ij/sqrt(self_similarity(i)*self_similarity(j));
            K(i,j) = k_ij;
            %K(j,i) = k_ij;
        end
    end
end

for i = 1:s
    K(i,i) = 1;
    for j = i+1:s
        K(j,i) = K(i,j);
    end
end
%{
parfor i = 1:s
    for j = 1:s
        k_ij = radial_structural_similarity(Ss{i}, Ss{j}, LSs{i}, LSs{j});
        k_ij = k_ij/sqrt(self_similarity(i)*self_similarity(j));
        K(i,j) = k_ij;
    end
end
%}

if verbose; disp('solving equation for alphas'); end
K = K.^zeta + lambda*eye(s);
alpha = K\energy;


%Prediction
if verbose; disp('making prediction'); end

if verbose; disp('computing spectra'); end
s_p = size(hold_out_set,1);
Ss_p = cell(s_p,1);
energy_p = zeros(s_p,1);
parfor idx = 1:s_p
    %disp(idx);
    h = hold_out_set(idx);
    N_p = molecule2neighbourhoods(h);
    Ss_p{idx} = neighbourhoods2spectra(N_p);
    energy_p(idx) = h.energy;
end

if verbose; disp('computing local similarity'); end
LSs_p = cell(s_p,1);    %local similarities
parfor idx = 1:s_p
    spectra = Ss_p{idx};
    n = size(spectra,1);
    LS = zeros(n,1);
    for i = 1:n
        LS(i) = radial_local_similarity(spectra(i,:),spectra(i,:));
    end
    LSs_p{idx} = LS;
end

if verbose; disp('computing self structural similarity'); end
self_similarity_p = zeros(s_p,1);  %for efficient normalisation
parfor idx = 1:s_p
    spectra_p = Ss_p{idx};
    LS = LSs_p{idx};
    self_similarity_p(idx) = radial_structural_similarity(spectra_p,spectra_p,LS,LS);
end

if verbose; disp('computing cross structural similarity'); end
L = zeros(s,s_p);
parfor i = 1:s
    for j = 1:s_p
        l_ij = radial_structural_similarity(Ss{i}, Ss_p{j}, LSs{i}, LSs_p{j});
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
