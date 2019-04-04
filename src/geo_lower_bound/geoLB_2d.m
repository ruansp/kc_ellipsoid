function [c_space, volPolyFit] = geoLB_2d(a, infla, th_b, N, isplot)
% geoLB_2d constructs KC space as a union of closed-form fitted polygon
% at different angles (Geometric Lower Bound)
%
% Inputs:
%  a     : semi-axis length of the smaller ellipse
%  infla : inflation factor ( >0 )
%  th_b  : rotational angle of the larger ellipse
%  N     : number of sampled points in a sphere
%  isplot: whether to plot the c-space
% 
% Outputs:
%  c_space   : vertices of configurations in C-space
%  volPolyFit: Volume of the fitted polyhedron
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

%% Parameters
alpha = a(1)/a(2);
b = a*(1+infla);

%% Max angle of rotation
phi = maxAngle(alpha,infla);

%% For each angle, compute mink diff
dphi = 0.01;
ang = -phi:dphi:phi;
c_space = []; area_poly = []; ang2 = [];
for i = 1:size(ang,2)
    % Closed-form Ellipse Fit at Shurnk Space
    aVec = [a;ang(i)];
    bVec = [b;th_b];
    [X_fit_poly, area] = cvxShpInMink( aVec, bVec, N);
    c_space = [c_space, [X_fit_poly; ang(i)*ones(1,size(X_fit_poly,2))]];
    
    if ~isempty(area)
        area_poly = [area_poly, area];
        ang2 = [ang2, ang(i)];
    end
    
end

volPolyFit = trapz(ang2, area_poly);
if isplot
    plot3(c_space(1,:),c_space(2,:),c_space(3,:),'r--','LineWidth',1.25)
end