
function descriptor = bag_of_bonds(molecule,max_size)

    %construct BoBs
    elements = {'H','C','N','O','S'};
    elements = sort(elements)';
    elems_no = size(elements,1);
    bonds = cell(0);
    for i = 1:elems_no
        for j = i:elems_no
            bond = strcat(elements{i},elements{j});
            bonds{size(bonds,2)+1} = bond;   %append
        end
    end
    
    bags_ondiag = repmat({Bag(max_size)},size(elements));
    bob_ondiag = containers.Map(elements,bags_ondiag);
    offdiag_max = (max_size-1)*max_size/2;
    bags_offdiag = repmat({Bag(offdiag_max)},size(bonds));
    bob_offdiag = containers.Map(bonds,bags_offdiag);
    
    %fill in the bags
    s = size(molecule.atoms.ff_coord,1);
    for i = 1:s
        at_i = char(molecule.atoms.types(i));
        z_i = get_atomic_number(at_i);
        R_i = angstrom2au(molecule.atoms.ff_coord(i,:));
        for j = 1:i
            if i == j
                val = 0.5*z_i^2.4;
                bob_ondiag(at_i) = bob_ondiag(at_i).add(val);
            else
                at_j = char(molecule.atoms.types(j));
                z_j = get_atomic_number(at_j);
                R_j = angstrom2au(molecule.atoms.ff_coord(j,:));
                val = z_i*z_j/norm(R_i-R_j);
                bond = sort(strcat(at_i,at_j));
                bob_offdiag(bond) = bob_offdiag(bond).add(val);
            end
        end
    end
    
    %sort each bag and rearange them into a column vector
    bonds_no = size(bonds,1);
    descriptor = zeros(bonds_no*offdiag_max + elems_no*max_size,1);
    last_desc = 0;
    
    for i = 1:size(elements,1)
        bag = bob_ondiag(elements{i});
        descriptor(last_desc+1:last_desc+max_size) = sort(bag.elements,'descend');
        last_desc = last_desc+max_size;
    end
    
    for i = 1:size(bonds,1)
        bag = bob_offdiag(bonds{i});
        descriptor(last_desc+1:last_desc+offdiag_max) = sort(bag.elements,'descend');
        last_desc = last_desc+offdiag_max;
    end
end

%element is a string
function z = get_atomic_number(element)
    elements = {'H','C','N','O','S'};
    vals = [1,6,7,8,16];
    atomic_numbers = containers.Map(elements,vals);
    z = atomic_numbers(element);
end

function au = angstrom2au(angstrom)
    au = angstrom/0.529;
end
