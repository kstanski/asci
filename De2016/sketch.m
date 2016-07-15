verbose = false;

if verbose; disp('training'); end
if verbose; disp('computing spectra'); end
s = size(training_set_proper,1);
Ss = cell(s,1);    %spectras
energy = zeros(s,1);
idx = 1;
disp(idx);
t = training_set_proper(idx);
N = molecule2neighbourhoods(t);
Ss{idx} = neighbourhoods2spectra(N);
energy(idx) = t.energy;

%Training
if verbose; disp('computing self similarity'); end
self_similarity = zeros(s,1);  %for efficient normalisation
idx = 1;
spectra = Ss{idx};
radial_structural_similarity(spectra,spectra)
