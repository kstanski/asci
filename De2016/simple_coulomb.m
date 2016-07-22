
function Coulomb = simple_coulomb(nhood,at_energy,max_size)
    
    s = size(nhood,1);
    %{
    C = zeros(s);
    for i = 1:s
        R_i = nhood(i,:);
        for j = 1:i
            if i == j
                C(i,i) = realmin;
            else
                R_j = nhood(j,:);
                m_ij = norm(R_i-R_j);
                C(i,j) = m_ij;
                C(j,i) = m_ij;
            end
        end
    end
    %}    

    %sort by row norm
    [~,idx] = sort(row_norm(nhood),'descend');
    C_sort = zeros(s,3);
    for i = 1:s
        %{
        for j = 1:s
            C_sort(i,j) = C(idx(i),idx(j));
        end
        %}
        C_sort(i,:) = nhood(idx(i),:);
    end
    %{
    %rearange and keep the lower triangular part only
    Coulomb = zeros(max_size*(max_size+1)/2, 1);
    coulomb_last = 0;
    for i = 1:s
        for j = 1:i
            coulomb_last = coulomb_last+1;
            Coulomb(coulomb_last) = C_sort(i,j);
        end
    end
    %}
    lack = max_size - s;
    Coulomb = padarray(C_sort,[lack 0],'post');
    Coulomb = Coulomb(:);
    if norm(Coulomb) ~= 0
        %Coulomb = Coulomb/norm(Coulomb);
    end
end

function N = row_norm(A)
    l = size(A,1);
    N = zeros(l,1);
    for i = 1:l
        N(i) = norm(A(i,:));
    end
end
