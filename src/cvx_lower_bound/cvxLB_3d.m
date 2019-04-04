function [Z_extreme, volPoly] = cvxLB_3d(a, infla, Hhc_3D)
% cvxLB_3d constructs polyhedron in C-Space for KC (Convex Lower Bound)
%
% Inputs:
%  a     : semi-axis length of the smaller ellipse
%  infla : inflation factor ( >0 )
%  Hhc_3D: symbolic H, h, c matrices
%
% Outputs:
%  Z_extreme: Extreme vertices that defines the convex polyhedron
%  volPoly  : Volume of the lower bound
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

%% Parameters
b = a * (1+infla);
A = diag(a.^(-2));
VA = diag(a);
B = diag(b.^(-2));

%% Extreme points with largest magnitude
disp('Extreme points with largest magnitude via convex optimization')
[Z_max, ~] = KC_Extreme_3d(a, infla, Hhc_3D);
Z_extreme = Z_max;

%% Extreme points at each axis
disp('Extreme points at each axis')
% Finding c, closed-form solution (maximum rotational angle)
ratio(3) = a(1)/a(2); ratio(1) = a(2)/a(3); ratio(2) = a(1)/a(3);
for i = 1:length(ratio)
    % Rotational axis
    axis_extreme(i) = maxAngle(ratio(i), infla);
    
    % Translational axis
    axis_extreme(i+length(ratio)) = b(i)-a(i);
end

Z_end = zeros(2*length(axis_extreme), length(axis_extreme));
for i = 1:length(axis_extreme)
    Z_end(2*i-1,i) = axis_extreme(i);
    Z_end(2*i,i) = -axis_extreme(i);
end

Z_extreme = [Z_extreme; Z_end];
%% Volume of polyhedron
volPoly = vol_poly(Z_extreme');