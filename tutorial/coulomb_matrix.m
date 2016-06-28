
function Coulomb = coulomb_matrix(molecule,max_size)
    s = size(molecule.atoms.ff_coord,1);
    C = zeros(s);
    for i = 1:s
        z_i = get_atomic_number(char(molecule.atoms.types(i)));
        R_i = angstrom2au(molecule.atoms.ff_coord(i,:));
        for j = 1:i
            if i == j
                C(i,i) = 0.5*z_i^2.4;
            else
                z_j = get_atomic_number(char(molecule.atoms.types(j)));
                R_j = angstrom2au(molecule.atoms.ff_coord(j,:));
                m_ij = z_i*z_j/norm(R_i-R_j);
                C(i,j) = m_ij;
                C(j,i) = m_ij;
            end
        end
    end
    
    %sort by row norm
    [~,idx] = sort(row_norm(C),'descend');
    C_sort = zeros(s);
    for i = 1:s
        for j = 1:s
            C_sort(i,j) = C(idx(i),idx(j));
        end
    end
    
    %rearange and keep the lower triangular part only
    Coulomb = zeros(max_size*(max_size+1)/2, 1);
    coulomb_last = 0;
    for i = 1:s
        for j = 1:i
            coulomb_last = coulomb_last+1;
            Coulomb(coulomb_last) = C_sort(i,j);
        end
    end
end

%element is a string
function z = get_atomic_number(element)
    elements = {'H','C','N','O','S'};
    vals = [1,6,7,8,16];
    atomic_numbers = containers.Map(elements,vals);
    z = atomic_numbers(element);
end

function N = row_norm(A)
    l = size(A,1);
    N = zeros(l,1);
    for i = 1:l
        N(i) = norm(A(i,:));
    end
end

function au = angstrom2au(angstrom)
    au = angstrom/0.529;
end
