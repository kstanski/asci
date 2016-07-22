
function S = neighbourhoods2spectra(N)
types_no = 5;    %H,C,N,O,S
S = cell(size(N,1),types_no);
for center_idx = 1:size(N,1)
    center = N{center_idx};
    for type = 1:types_no    
        S{center_idx,type} = power_spectrum(center{type});
        %S{center_idx,type} = simple_coulomb(center{type},get_at_energy(type),16);
    end
end
end

function at_energy = get_at_energy(idx)
Z = [1,6,7,8,16];
z_i = Z(idx);
at_energy = 0.5*z_i^2.4;
end
