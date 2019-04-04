% Comparisons of C-space volume between different methods in 3D
% 
% Variables:
%   Aspect Ratios, i.e. a1/a2, a1/a3
% 
% C-space to be compared:
%  (1) Convex Lower Bound
%  (2) Geometric Lower Bound
%  (3) "Actual" KC C-space: from Minkowski difference
%
% Plots:
%  (1) Absolute volumes in C-space
%  (2) Relative volumes in C-space with respect to the "Actual" KC C-space
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

%% Initialization
clc; clear; close all;

% add paths
addpath ../include/
addpath ../mat/
addpath ../src/cvx_lower_bound/
addpath ../src/geo_lower_bound/

load('Hhc_3D.mat')
infla = 0.1;
a(1,1) = 4;

Nasp = 5;
amax = 1.2;
amin = 1.001;
[asp1,asp2] = meshgrid(amin:(amax-amin)/(Nasp-1):amax,...
    amin:(amax-amin)/(Nasp-1):amax);

%% Construct c-space based on different methods
% figure; hold on; axis equal; grid on;
iter = 0;
for i = 1:size(asp1,1)
    for j = 1:size(asp1,2)
        iter = iter+1;
        disp(['Loop: ', num2str(iter), ', Aspect Ratio: ',...
            num2str(asp1(i,j)), ', ', num2str(asp2(i,j))]);
        
        a(2,1) = a(1,1)/asp1(i,j);
        a(3,1) = a(1,1)/asp2(i,j);
        
        b = a*(1+infla);
        % Convex Lower Bound
        disp('==== Convex Lower Bound ====')
        [Z_extreme, volPoly] = cvxLB_3d(a, infla, Hhc_3D);
        
        % Geometric Lower Bound
        disp('==== Geometric Lower Bound ====')
        [c_polyFit, volPolyFit] = geoLB_3d(a, infla, 10);
        
        % "Actual" KC C-space
        disp('==== "Actual" KC C-space ====')
        [c_space, volMink] = mink_3d(a, infla, 10);
        
        vp(i,j) = volPoly;
        vpf(i,j) = volPolyFit;
        vm(i,j) = volMink;
    end
end

%% Plots
figure; hold on; grid on;
fa = 0.3;
ea = 0.1;
surf(asp1, asp2, vm, 'FaceColor','k','FaceAlpha',fa,'EdgeAlpha',ea);
surf(asp1, asp2, vp, 'FaceColor','b','FaceAlpha',fa,'EdgeAlpha',ea);
surf(asp1, asp2, vpf, 'FaceColor','r','FaceAlpha',fa,'EdgeAlpha',ea);
legend('"Actual" KC C-space', 'Convex Lower Bound', 'Geometric Lower Bound')
xlabel('\alpha_1')
ylabel('\alpha_2')
zlabel('Volume')

figure; hold on; grid on;
fa = 0.3;
ea = 0.1;
surf(asp1, asp2, vm./vm, 'FaceColor','k','FaceAlpha',fa,'EdgeAlpha',ea);
surf(asp1, asp2, vp./vm, 'FaceColor','b','FaceAlpha',fa,'EdgeAlpha',ea);
surf(asp1, asp2, vpf./vm, 'FaceColor','r','FaceAlpha',fa,'EdgeAlpha',ea);
legend('"Actual" KC C-space', 'Convex Lower Bound', 'Geometric Lower Bound')
xlabel('\alpha_1')
ylabel('\alpha_2')
zlabel('Relative Volume')

