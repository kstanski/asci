
tsp = training_set_proper(1);
n_at = size(tsp.atoms.ff_coord,1);
%permutation
%tsp.atoms.ff_coord = vertcat(tsp.atoms.ff_coord(3:n_at,:),tsp.atoms.ff_coord(1:2,:));
%tsp.atoms.types = vertcat(tsp.atoms.types(3:n_at,:),tsp.atoms.types(1:2,:));
%translation
%tsp.atoms.ff_coord = tsp.atoms.ff_coord + 1.5;
%rotation
%tsp.atoms.ff_coord = quatrotate([pi/0.9,pi/2.2,pi/3.01,pi/4.1],tsp.atoms.ff_coord);
for i = 1:n_at
    x = tsp.atoms.ff_coord(i,1);
    y = tsp.atoms.ff_coord(i,2);
    z = tsp.atoms.ff_coord(i,3);
    [a,b,c] = cart2sph(x,y,z);
    a = a+0.56;
    b = b+1.63;
    [x,y,z] = sph2cart(a,b,c);
    tsp.atoms.ff_coord(i,:) = [x,y,z];
end

%arbitrary alteration
tsp.atoms.ff_coord(2,3) = tsp.atoms.ff_coord(2,3) + 1;

N1 = molecule2neighbourhoods(training_set_proper(1));
N2 = molecule2neighbourhoods(tsp);

S1 = neighbourhoods2spectra(N1);
S2 = neighbourhoods2spectra(N2);

tic
ss = radial_structural_similarity(S1,S2);
ss = ss/sqrt(radial_structural_similarity(S1,S1)*radial_structural_similarity(S2,S2));
disp(sprintf('%.10f',ss))
toc