
function [RMSE,MAE,R2] = basic_model(lambda,sigma,training_set_proper,hold_out_set,kernel)

%Training
n = length(training_set_proper);
X = cell(n,1);
Y = zeros(n,1);

for i = 1:n
    T = training_set_proper(i);
    X{i} = compute_coulomb_matrix(T,23);
    Y(i) = T.energy;
end

%lambda = 10^(-6.5);
%sigma = 724;

K = zeros(n);
for i = 1:n
    for j = 1:i
        k_ij = kernel(X{i}, X{j}, sigma);
        K(i,j) = k_ij;
        K(j,i) = k_ij;
    end
end

K = K + lambda*eye(n);
%U = chol(K);    %Cholesky decomposition
%alpha = forward_backward_substitution(U,Y);
alpha = K\Y;


%Prediction
m_p = length(hold_out_set);
X_p = cell(m_p,1);
Y_p = zeros(m_p,1);

for i = 1:m_p
    T_p = hold_out_set(i);
    X_p{i} = compute_coulomb_matrix(T_p,23);
    Y_p(i) = T_p.energy;
end

L = zeros(n,m_p);
for i = 1:n
    for j = 1:m_p
        L(i,j) = kernel(X{i}, X_p{j}, sigma);
    end
end

f = L.'*alpha;

RMSE = rms(Y_p-f);
MAE = mean(abs(Y_p-f));
R2 = corr(Y_p,f)^2;

a = linspace(-2200,-1000,2);
hold on
scatter(Y_p,f,10,'filled');
plot(a,a,'k');

end

