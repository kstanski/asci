
function X = rotate_invariant(X)
s = size(X,1);

%move to 'mass' center
%R_c = mean(X);
%X = X - ones(s,1)*R_c;

%{
%rotate by average angle towards positive z axis (phi)
[~,Phi,~] = spherical_coords(X);
q_sum = zeros(1,4);    %summative quaternion
for idx = 1:s
    z = [0,0,1];
    v = X(idx,:);
    r_axis = cross(z,v);    %direction of rotation axis
    if norm(r_axis) ~= 0
        r_axis = r_axis/norm(r_axis);
        q = get_quat(Phi(idx),r_axis);
        q_sum = q_sum + q;
    end
end
q_avg = quatnormalize(q_sum);
for idx = 1:s
    X(idx,:) = quatrotate(q_avg,X(idx,:));
end
%}
%{
%rotate by average angle towards zy plane (theta)
[Theta,~,~] = spherical_coords(X);
t_m = mean(Theta);
q_sum = zeros(1,4);    %summative quaternion
for idx = 1:s
    z = [0,0,1];
    q = get_quat(Theta(idx),z);
    q_sum = q_sum + quatnormalize(q);
end
q_avg = quatnormalize(q_sum);
%q_avg = get_quat(t_m,[0,0,1]);
for idx = 1:s
    X(idx,:) = quatrotate(q_avg,X(idx,:));
end
[Theta,~,~] = spherical_coords(X);
t_m = mean(abs(Theta));

X_180 = zeros(s,3);
for idx = 1:s
    X_180(idx,:) = quatrotate(get_quat(pi,[0,0,1]),X(idx,:));
end
[Theta,~,~] = spherical_coords(X_180);
t_m_180 = mean(abs(Theta));

if t_m_180 < t_m
    X = X_180;
end
%}
if size(X,2) == 3
    %allign with z axis (phi)
    X_mod = zeros(size(X));
    for idx = 1:s
        x_i = X(idx,:);
        X_mod(idx,:) = x_i*norm(x_i);
    end
    R_c = mean(X_mod,1);
    if R_c == 0
        if X(1,:) ~= [0,0,0]
            disp('not center')
        end
        v = X(1,:);
    else
        v = R_c;
    end
    z = [0,0,1];
    axis = cross(v,z);
    if norm(axis) ~= 0
        axis = axis/norm(axis);
    end
    angle = get_angle(v,z);
    q_allign_x = get_quat(-angle,axis);
    for idx = 1:s
        X(idx,:) = quatrotate(q_allign_x,X(idx,:));
    end
    
    %allign with zx plane (theta)
    R_c = mean(X,1);
    v = [R_c(1),R_c(2),0];
    x = [1,0,0];
    axis = cross(v,x);
    if norm(axis) ~= 0
        axis = axis/norm(axis);
    end
    angle = get_angle(v,x);
    q_allign_x = get_quat(-angle,axis);
    for idx = 1:s
        X(idx,:) = quatrotate(q_allign_x,X(idx,:));
    end
end 
    
end

function q = get_quat(alpha,R)
s = sin(alpha/2);
q = [cos(alpha/2),s*R(1),s*R(2),s*R(3)];
end

function angle = get_angle(a,b)
angle = atan2(norm(cross(a,b)),dot(a,b));
end

function [Theta,Phi,R] = spherical_coords(X)
if size(X,2) == 3
    [Theta,Elev,R] = cart2sph(X(:,1),X(:,2),X(:,3));
else
    Theta = 0;
    Elev = 0;
    R = 0;
end

Phi = 2*pi - Elev;
end
