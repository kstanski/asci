
lambda = 2.^[-40:0.5:-5];
sigma = 2.^[5:0.5:18];

grid = cell(length(lambda),length(sigma));

for l = 1:length(lambda)
    for s = 1:length(sigma)
        a = zeros(3,1);
        [a(1),a(2),a(3)] = basic_model(lambda(l),sigma(s),training_set_proper,hold_out_set);
        grid{l,s} = a;
    end
end