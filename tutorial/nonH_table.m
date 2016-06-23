
%molecules = read_data('dsgdb7ae2.xyz',7102);
count_nonH = zeros(8,1);

for i = 1:length(molecules)
    m = molecules(i);
    c = m.atoms.nonH;
    if c > 0
        count_nonH(c) = count_nonH(c)+1;
    end
end

count_nonH(8) = sum(count_nonH(1:7));
