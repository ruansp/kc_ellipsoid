function [Z_extreme, volPoly] = cvxLB_2d(a, infla, th_b, Hhc_rot_2D, isplot)
% cvxLB_2d constructs KC space as a polyhedron based on extreme vertices 
%  (Convex Lower Bound)
%
% Inputs:
%  a     : semi-axis length of the smaller ellipse
%  infla : inflation factor ( >0 )
%  th_b  : rotational angle of the larger ellipse
%  Hhc_2D: symbolic H, h, c matrices
%  isplot: whether to plot the c-space
%
% Outputs:
%  Z_extreme: Extreme vertices that defines the convex polyhedron
%  volPoly  : Volume of the lower bound
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019


%% Extreme points with farthest distance to origin
Z_max = [];
while isempty(Z_max)
    [Z_max, ~] = KC_Extreme(a, infla, th_b, Hhc_rot_2D);
end

Z_extreme(:,1) = Z_max(:,2);
Z_extreme(:,2) = Z_max(:,3);
Z_extreme(:,3) = Z_max(:,1);

x = Z_extreme(:,1); y = Z_extreme(:,2); z = Z_extreme(:,3);

%% Extreme points at each axis
% Finding c, closed-form solution (maximum rotational angle)
ratio = a(1)/a(2);
c = maxAngle(ratio, infla);

% Extreme value at each axis
aa = a(1)*infla;
bb = a(2)*infla;
cc = c;

Z_end = [aa,0,0; -aa,0,0;...
    0,bb,0; 0,-bb,0;...
    0,0,cc; 0,0,-cc];

Z_extreme = [Z_extreme; Z_end];
x = Z_extreme(:,1); y = Z_extreme(:,2); z = Z_extreme(:,3);

%% Volume of the polyhedron
volPoly = vol_poly(Z_extreme');
K = convhull(x,y,z);

%% Plots
if isplot
    plot3(x,y,z,'b*', 'LineWidth', 2);
    patch('Faces',K, 'Vertices',[x y z], 'FaceColor', 'y', 'FaceAlpha', 0.3);
end