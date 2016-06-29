
lambda = 2.^[-40:2:-30];
sigma = 2.^[35:0.5:45];

kernel = @laplacian_kernel;
descriptor = @bag_of_bonds;
verbose = false;

[X,Y,X_p,Y_p] = apply_descriptor(training_set_proper,hold_out_set,descriptor);

grid = repmat(100,length(lambda),length(sigma),3);
for l = 1:length(lambda)
    for s = 1:length(sigma)
        disp(sprintf('lambda = %e\nsigma = %f\n',lambda(l),sigma(s)))
        a = zeros(3,1);
        [~,a(1),a(2),a(3)] = krr(lambda(l),sigma(s),X,Y,X_p,Y_p,kernel,verbose);
        grid(l,s,:) = a;
    end
end

mae_arr = grid(:,:,2);
[M,I] = min(mae_arr(:));
[min_l,min_s] = ind2sub(size(mae_arr),I);

disp(sprintf('\nmin MAE = %f for:\nlambda = %e\nsigma = %f',mae_arr(min_l,min_s),lambda(min_l),sigma(min_s)))
