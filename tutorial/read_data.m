
function molecules = read_data(filename,N)
    fid = fopen(filename,'r');

    molecules = repmat(Molecule(0), N, 1);

    count_mol = 0;
    count_at = 0;
    line = fgetl(fid);
    while ischar(line)
        splt_line = strsplit(line);
        switch length(splt_line)
            case 0
            case 1
                count_mol = count_mol+1;
                count_at = 0;
                atoms_no = str2double(splt_line(1));
                molecule = Molecule(atoms_no);
            case 2
                molecule.id = splt_line(1);
                molecule.energy = str2double(splt_line(2));
            otherwise
                count_at = count_at+1;
                at = splt_line(1);
                if ~strcmp(at,'H')
                    molecule.size = molecule.size + 1;
                end
                molecule.atoms.types(count_at) = at;
                molecule.atoms.ff_coord(count_at,:) = str2double(splt_line(2:4));
                molecule.atoms.dft_coord(count_at,:) = str2double(splt_line(5:7));
                if eq(count_at,atoms_no)
                    molecules(count_mol) = molecule;
                end      
        end
        line = fgetl(fid);
    end

    fclose(fid);
end
