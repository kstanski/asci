
r_cut = 3;
R = linspace(0,r_cut,100);
G = zeros(length(R),1);

hold on
for n = 1:4
    for i = 1:length(R)
        G(i) = radial_basis_function(R(i),r_cut,n);
    end
    plot(R,G);
end
hold off

