
train_no = 1000;
all_no = 7102;

molecules = read_data('dsgdb7ae2.xyz',all_no);
%molecules = read_data('test_data.xyz',2);
below5atoms = repmat(Molecule(0), train_no, 1);
below5_last = 0;
above5atoms = repmat(Molecule(0), all_no, 1);
above5_last = 0;

for i = 1:length(molecules)
    m = molecules(i);
    if m.size < 5
        below5_last = below5_last+1;
        below5atoms(below5_last) = m;
    else
        above5_last = above5_last+1;
        above5atoms(above5_last) = m;
    end
end

below5atoms = below5atoms(1:below5_last);
%sort the set before stratification
above5atoms = above5atoms(1:above5_last);
[~,idx] = sort([above5atoms.size]);
above5atoms = above5atoms(idx);

prediction_set = repmat(Molecule(0), all_no-train_no, 1);
training_set = repmat(Molecule(0), train_no-below5_last, 1);    %for above5 only!
[training_set,prediction_set] = stratify(above5atoms,training_set,prediction_set);

%Create a hold-out set
hold_out_set = repmat(Molecule(0), 100, 1);
training_set_proper = repmat(Molecule(0), 900-length(below5atoms), 1);

[hold_out_set,training_set_proper] = stratify(training_set,hold_out_set,training_set_proper);
training_set_proper = cat(1,below5atoms,training_set_proper);
