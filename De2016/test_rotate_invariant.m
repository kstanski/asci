N1 = molecule2neighbourhoods(training_set_proper(50));
n1 = N1{1};

X = n1{1};

%X = [12,65,34;1,2,1.85;23,63,47]

rx = rotate_invariant(X);

hold on
scatter(rx(:,1),rx(:,2),'bx')
%scatter(X(:,1),X(:,2),'b.')
X = quatrotate([0.5,2,3,4],X);

n1 = X;
for i = 1:size(n1,1)
    x = n1(i,1);
    y = n1(i,2);
    z = n1(i,3);
    [a,b,c] = cart2sph(x,y,z);
    a = a+pi/3.56;
    b = b+0;
    c = c+0;
    [x,y,z] = sph2cart(a,b,c);
    n1(i,:) = [x,y,z];
end

rn = rotate_invariant(n1);
scatter(rn(:,1),rn(:,2),'g')
%scatter(n1(:,1),n1(:,2),'g.')

rn - rx