function [c_space, volMink] = mink_2d(a, infla, th_b, N, isplot)
% mink_2d constructs KC space as a union of Minkowski difference at
% different angles
%
% Inputs:
%  a     : semi-axis length of the smaller ellipse
%  infla : inflation factor ( >0 )
%  th_b  : rotational angle of the larger ellipse
%  N     : number of points on the Minkowski difference boundary
%  isplot: whether to plot the c-space
%
% Dependency:
%  "Ellipsoidal Toolbox": version 1.2.0.0 by Alex Kurzhanskiy
%  URL: https://www.mathworks.com/matlabcentral/fileexchange/21936-ellipsoidal-toolbox-et
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

%% Parameters
alpha = a(1)/a(2);
b = a*(1+infla);

% Initialize ET parameters
global ellOptions;
ellOptions.verbose = 0;
ellOptions.rel_tol = 1e-6;
ellOptions.abs_tol = 1e-7;
ellOptions.plot2d_grid = N;

VA = diag(a.^2);
VB = diag(b.^2);

%% Max angle of rotation
phi = maxAngle(alpha,infla);

%% For each angle, compute mink diff
ang = -phi:0.01:phi;
c_space = []; area_mink = []; ang2 = [];
for i = 1:size(ang,2)
    % Closed-form Mink difference
    aVec = [a;ang(i)];
    bVec = [b;th_b];
    
    % Ellipsoid using ET
    A = rot2(ang(i)) * VA * rot2(ang(i))';
    B = rot2(th_b) * VB * rot2(th_b)';
    
    A(1,2) = A(2,1);
    B(1,2) = B(2,1);
    
    Ea = ellipsoid(A);
    Eb = ellipsoid(B);
    
    % Compute Minkowski difference
    [~, X_eb] = minkdiff(Eb, Ea);
    c_space = [c_space, [X_eb; ang(i)*ones(1,size(X_eb,2))]];
    
    if size(X_eb,2) >= 3
        [~,area] = convhull(X_eb');
        area_mink = [area_mink, area];
        ang2 = [ang2, ang(i)];
    end
end

volMink = trapz(ang2, area_mink);
if isplot
    plot3(c_space(1,:),c_space(2,:),c_space(3,:),'k.');
end