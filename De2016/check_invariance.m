
tsp = training_set_proper(2);
n_at = size(tsp.atoms.ff_coord,1);
%permutation
%tsp.atoms.ff_coord = vertcat(tsp.atoms.ff_coord(3:n_at,:),tsp.atoms.ff_coord(1:2,:));
%tsp.atoms.types = vertcat(tsp.atoms.types(3:n_at,:),tsp.atoms.types(1:2,:));
%translation
%tsp.atoms.ff_coord = tsp.atoms.ff_coord + 1.5;
%rotation
%tsp.atoms.ff_coord = quatrotate([pi/0.9,pi/2.2,pi/3.01,pi/4.1],tsp.atoms.ff_coord);
%{
for i = 1:n_at
    x = tsp.atoms.ff_coord(i,1);
    y = tsp.atoms.ff_coord(i,2);
    z = tsp.atoms.ff_coord(i,3);
    [a,b,c] = cart2sph(x,y,z);
    a = a+0;
    b = b+0;
    c = c+0;
    [x,y,z] = sph2cart(a,b,c);
    tsp.atoms.ff_coord(i,:) = [x,y,z];
end
%}

%arbitrary alteration
%tsp.atoms.ff_coord(2,3) = tsp.atoms.ff_coord(2,3) + 1;

N1 = molecule2neighbourhoods(training_set_proper(2));
N2 = molecule2neighbourhoods(tsp);

n1 = N1{1};
n1 = n1{1};
n1 = quatrotate([pi/0.9,pi,pi/3.01,pi/4.1],n1);

for i = 1:size(n1,1)
    x = n1(i,1);
    y = n1(i,2);
    z = n1(i,3);
    [a,b,c] = cart2sph(x,y,z);
    a = a+0;
    b = b+0;
    c = c+0;
    [x,y,z] = sph2cart(a,b,c);
    n1(i,:) = [x,y,z];
end

%n1 = rotate_invariant(n1);
ps1 = power_spectrum(n1);
n2 = N2{1};
%n2{1} = rotate_invariant(n2{1});
ps2 = power_spectrum(n2{1});
%soap(n1,n2{1})/sqrt(soap(n1,n1)*soap(n2{1},n2{1}))
dot(ps1,ps2)

%S1 = neighbourhoods2spectra(N1);
%S2 = neighbourhoods2spectra(N2);

%tic
%ss = radial_structural_similarity(S1,S2);
%ss = ss/sqrt(radial_structural_similarity(S1,S1)*radial_structural_similarity(S2,S2));
%disp(sprintf('%.10f',ss))
%toc