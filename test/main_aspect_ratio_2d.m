% Comparisons of C-space volume between different methods in 2D
% 
% Variables:
%   Aspect Ratios, i.e. a1/a2
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

load('Hhc_rot_2D.mat')
infla = 0.08;
alpha = 1.1:0.01:1.5;

for i = 1:size(alpha,2)
% for i = 3
    disp(['Loop: ', num2str(i), ', Aspect Ratio: ', num2str(alpha(i))]);
    a(1,1) = 5;
    a(2,1) = a(1)/alpha(i);
    
    b = a*infla;
    E1 = diag(a.^(-2));
    E2 = diag(b.^(-2));
    
    %% Construct c-space based on different methods 
    % Convex Lower Bound
    disp('==== Convex Lower Bound ====')
    [Z_extreme, volPoly] = cvxLB_2d(a, infla, 0, Hhc_rot_2D, 0);
    
    % Geometric Lower Bound
    disp('==== Geometric Lower Bound ====')
    [c_space3, volPolyFit] = geoLB_2d(a, infla, 0, 1e3, 0);
    
    % "Actual" KC C-space
    disp('==== "Actual" KC C-space ====')
    [c_space, volMink] = mink_2d(a, infla, 0, 1e3, 0);

    vp(i) = volPoly;
    vpf(i) = volPolyFit;
    vm(i) = volMink;
end

figure; hold on; grid on;
lw = 1.25;
plot(alpha, vm, 'k-', 'LineWidth', lw);
plot(alpha, vp, 'b-.', 'LineWidth', lw);
plot(alpha, vpf, 'r-*', 'LineWidth', lw);
legend('"Actual" KC C-space', 'Convex Lower Bound', 'Geometric Lower Bound')
xlabel('Aspect Ratio')
ylabel('Volume')

figure; hold on; grid on;
lw = 1.25;
plot(alpha, vm./vm, 'k-', 'LineWidth', lw);
plot(alpha, vp./vm, 'b-.', 'LineWidth', lw);
plot(alpha, vpf./vm, 'r-*', 'LineWidth', lw);
legend('"Actual" KC C-space', 'Convex Lower Bound', 'Geometric Lower Bound')
xlabel('Aspect Ratio')
ylabel('Relative Volume')
