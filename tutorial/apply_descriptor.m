function [X,Y,X_p,Y_p] = apply_descriptor(training_set_proper,hold_out_set,descriptor)

disp('computing training descriptor vectors')
n = length(training_set_proper);
X = cell(n,1);
Y = zeros(n,1);

for i = 1:n
    T = training_set_proper(i);
    X{i} = descriptor(T,23);
    Y(i) = T.energy;
end

disp('computing evaluation descriptor vectors')
m_p = length(hold_out_set);
X_p = cell(m_p,1);
Y_p = zeros(m_p,1);

for i = 1:m_p
    T_p = hold_out_set(i);
    X_p{i} = descriptor(T_p,23);
    Y_p(i) = T_p.energy;
end

end
