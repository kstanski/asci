
function S = neighbourhoods2spectra(N)
types_no = 5;    %H,C,N,O,S
S = cell(size(N,1),types_no);
for center_idx = 1:size(N,1)
    center = N{center_idx};
    for type = 1:types_no    
        S{center_idx,type} = power_spectrum(center{type});
    end
end
end
