classdef Molecule    
    properties
        id
        energy
        atoms
        size    %number of non-H atoms
    end
    
    methods
        function obj = Molecule(atoms_no)
           obj.atoms.types = repmat({' '}, atoms_no, 1);
           obj.atoms.ff_coord = zeros(atoms_no,3);
           obj.atoms.dft_coord = zeros(atoms_no,3);
           obj.size = 0;
        end
    end
    
end

