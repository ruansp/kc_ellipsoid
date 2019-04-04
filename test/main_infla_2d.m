% Comparisons of C-space volume between different methods in 2D
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

%% Initialization
clc; clear; close all;

% add paths
addpath ../include/
addpath ../mat/
addpath ../src/cvx_lower_bound/
addpath ../src/geo_lower_bound/

%% Parameters
load('Hhc_rot_2D.mat')
infla = 0.05:0.005:0.1;
a = [5;3.5];

%% Construct c-space based on different methods
for i = 1:size(infla,2)
    disp(['Loop: ', num2str(i), ', Inflation Factor: ', num2str(infla(i))]);
    
    % Convex Lower Bound
    disp('==== Convex Lower Bound ====')
    [Z_extreme, volCvx] = cvxLB_2d(a, infla(i), 0, Hhc_rot_2D, 0);
    
    % Geometric Lower Bound
    disp('==== Geometric Lower Bound ====')
    [c_space3, volGeo] = geoLB_2d(a, infla(i), 0, 1e3, 0);
    
    % Actual KC C-space: Minkowski difference
    disp('==== Actual KC C-space: Minkowski Difference ====')
    [c_space, volMink] = mink_2d(a, infla(i), 0, 1e3, 0);
    
    vp(i) = volCvx;
    vpf(i) = volGeo;
    vm(i) = volMink;
end

figure; hold on; grid on;
lw = 1.25;
plot(infla, vm, 'k-', 'LineWidth', lw);
plot(infla, vp, 'b-.', 'LineWidth', lw);
plot(infla, vpf, 'r--', 'LineWidth', lw);
legend('Actual KC C-space', 'Convex Lower Bound',...
    'Geometric Lower Bound')
xlabel('Inflation Factor')
ylabel('Volume')

figure; hold on; grid on;
lw = 1.25;
plot(infla, vm./vm, 'k-', 'LineWidth', lw);
plot(infla, vp./vm, 'b-.', 'LineWidth', lw);
plot(infla, vpf./vm, 'r--', 'LineWidth', lw);
legend('Actual KC C-space', 'Convex Lower Bound',...
    'Geometric Lower Bound')
xlabel('Inflation Factor')
ylabel('Relative Volume')

