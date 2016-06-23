
kernel = @laplacian_kernel;

[RMSE,MAE,R2] = basic_model(10^(-6.5),724,training_set_proper,hold_out_set,kernel);
disp(sprintf('RMSE: %f\nMAE: %f\nR2: %f\n',RMSE,MAE,R2))
