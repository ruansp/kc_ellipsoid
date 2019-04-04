% Comparisons of C-space volume between different methods in 3D
% 
% Variables:
%   Inflation Factor, i.e. infla = b/a-1
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

clc; clear; close all;

%% Parameters
addpath ../include/
addpath ../mat/
addpath ../src/cvx_lower_bound/
addpath ../src/geo_lower_bound/

load('Hhc_3D.mat')
infla = 0.005:0.005:0.1;
a = [4;2.5;2];

i = 1;
b = a*(1+infla(i));
E1 = diag(a.^(-2));
E2 = diag(b.^(-2));

%% Construct c-space based on different methods
for i = 1:size(infla,2)
    disp(['Loop: ', num2str(i), ', Inflation Factor: ', num2str(infla(i))]);
    
    % Convex Lower Bound
    disp('==== Convex Lower Bound ====')
    [Z_extreme, volPoly] = cvxLB_3d(a, infla(i), Hhc_3D);
    
    % Geometric Lower Bound
    disp('==== Geometric Lower Bound ====')
    [c_polyFit, volPolyFit] = geoLB_3d(a, infla(i), 10);
    
    % "Actual" KC C-space
    disp('==== "Actual" KC C-space ====')
    [c_space, volMink] = mink_3d(a, infla(i), 10);

    vp(i) = volPoly;
    vpf(i) = volPolyFit;
    vm(i) = volMink;
end

%% Plots
figure; hold on; grid on;
lw = 1.25;
plot(infla, vm, 'k-', 'LineWidth', lw);
plot(infla, vp, 'b-.', 'LineWidth', lw);
plot(infla, vpf, 'r--', 'LineWidth', lw);
legend('"Actual" KC C-space', 'Convex Lower Bound',...
    'Geometric Lower Bound')
xlabel('Inflation Factor')
ylabel('Volume')

figure; hold on; grid on;
lw = 1.25;
plot(infla, vm./vm, 'k-', 'LineWidth', lw);
plot(infla, vp./vm, 'b-.', 'LineWidth', lw);
plot(infla, vpf./vm, 'r--', 'LineWidth', lw);
legend('"Actual" KC C-space', 'Convex Lower Bound',...
    'Geometric Lower Bound')
xlabel('Inflation Factor')
ylabel('Relative Volume')

