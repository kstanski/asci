function [RMSE,MAE,R2] = plot_at_energy(Y_p,f)

RMSE = rms(Y_p-f);
MAE = mean(abs(Y_p-f));
R2 = corr(Y_p,f)^2;

a = linspace(-2200,-1000,2);
hold on
scatter(Y_p,f,10,'filled');
plot(a,a,'k');
xlabel('reference (kcal/mol)');
ylabel('predicted (kcal/mol)');
stats = sprintf('RMSE: %f\nMAE: %f\nR2: %f\n',RMSE,MAE,R2);
text(-1500,-1900,stats);
hold off

end