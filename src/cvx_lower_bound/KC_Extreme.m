function [Z_max, f_max] = KC_Extreme(a, infla, angle, Hhc_rot_2D)
% KC_Extreme computes the extreme vertices for Convex Lower Bound in 2D
%
% Inputs:
%  a         : semi-axis length of the smaller ellipse
%  infla     : inflation factor ( >0 )
%  angle     : Angle of the 2D rotation for the fixed ellipsoid
%  Hhc_rot_2D: symbolic H, h, c matrices
%
% Outputs:
%  Z_max: Extreme configuration with the maximum magnitude
%  f_max: The corresponding value of the cost function
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

infla = infla+1;
%% Initialization: numerical H, h, c
syms a1 a2 b1 b2 u1 u2 theta
syms th x y
value = [a(1), a(2), a(1)*infla, a(2)*infla, angle];

H = subs( Hhc_rot_2D.H, [a1, a2, b1, b2, theta], value);
h = subs( Hhc_rot_2D.h, [a1, a2, b1, b2, theta], value);
c = subs( Hhc_rot_2D.c, [a1, a2, b1, b2, theta], value);
H = vpa(H, 4);
h = vpa(h, 4);
c = vpa(c, 4);

n_start = 100; % Number of starting guessing points
dtheta = pi/6/n_start;
theta = 0:dtheta:pi/2;
z0 = [theta; rand(2, size(theta,2))];

%% Objective function
z = [th; x; y];
f = -z'*z;
f2 = matlabFunction(f,'vars',{z});

%% Inequality Constraint
ang = 0:pi/20:2*pi;
u = [cos(ang); sin(ang)];

for i = 1:size(ang,2)
    H2 = subs( H, [u1, u2], u(1:2,i)');
    h2 = subs( h, [u1, u2], u(1:2,i)');
    c2 = subs( c, [u1, u2], u(1:2,i)');
    H2 = double(H2);
    h2 = double(h2);
    c2 = double(c2);
    cieq(i) = z'*H2*z + h2'*z + c2 - 1;
end

constraint = matlabFunction(cieq, [], 'vars', {z});

%% Optimization
options = optimoptions('fmincon','Algorithm','sqp','Display','off');
for i = 1:length(theta)
    [z_x(:,i), fval(i)] = fmincon(f2,z0(:,i),[],[],[],[],[],[],constraint,options);
end

[f_max, idx] = max(-fval);
z_max = z_x(:,idx);

%% Validation of possible solutions
% 8 possible solutions
Z_sol = [ z_max(1), z_max(2),  z_max(3);...
    z_max(1), z_max(2), -z_max(3);...
    z_max(1),-z_max(2),  z_max(3);...
    z_max(1),-z_max(2), -z_max(3);...
    -z_max(1), z_max(2),  z_max(3);...
    -z_max(1), z_max(2), -z_max(3);...
    -z_max(1),-z_max(2),  z_max(3);...
    -z_max(1),-z_max(2), -z_max(3)];

% Validation
Z_max = [];
for i = 1:size(Z_sol,1)
    inside = 1;
    
    % condition that smaller ellipse moves inside larger one without
    % collision.
    A = [a(1) 0; 0 a(2)];
    B = [a(1)*infla 0; 0 a(2)*infla];
    IneqF = (rot2(Z_sol(i,1))*A*u + Z_sol(i,2:3)')' * B^(-2) * (rot2(Z_sol(i,1))*A*u + Z_sol(i,2:3)');
    ineqF = diag(IneqF);
    
    inside = all(ineqF <= 1);
    
    
    if inside
        Z_max = [Z_max; Z_sol(i,:)];
    end
end
end