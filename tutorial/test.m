
kernel = @test_kernel;

[RMSE,MAE,R2] = basic_model(10^(-6.5),724,training_set_proper,hold_out_set,kernel)
