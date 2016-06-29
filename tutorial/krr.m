
function [f,RMSE,MAE,R2] = krr(lambda,sigma,X,Y,X_p,Y_p,kernel,verbose)

%Training
n = size(X,1);
if verbose; disp('kernelising'); end
K = zeros(n);
for i = 1:n
    for j = 1:i
        k_ij = kernel(X{i}, X{j}, sigma);
        K(i,j) = k_ij;
        K(j,i) = k_ij;
    end
end

if verbose; disp('solving equation for alphas'); end
K = K + lambda*eye(n);
alpha = K\Y;


%Prediction
if verbose; disp('making prediction'); end
m_p = size(X_p,1);
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

if verbose; plot_at_energy(Y_p,f); end

end
