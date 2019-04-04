% Visualizations of C-space for different methods in 2D
%
% C-space to be compared:
%  (1) Convex Lower Bound
%  (2) Geometric Lower Bound
%  (3) "Actual" KC C-space: from Minkowski difference
%
% Plots:
%  (1) The three C-spaces in SE(2) configuration space
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

clc; clear; close all;

%% Parameters
addpath ../mat/
addpath ../src/cvx_lower_bound/
addpath ../src/geo_lower_bound/
addpath ../include/

load('Hhc_rot_2D.mat')
infla = 0.1;
alpha = 5/3.5;

a(1,1) = 5;
a(2,1) = a(1)/alpha;

b = a*infla;
E1 = diag(a.^(-2));
E2 = diag(b.^(-2));

%% Construct c-space based on different methods
figure; hold on; axis equal; grid on;
% Polyhegron
disp('==== Convex Lower Bound ====')
[Z_extreme, volPoly] = cvxLB_2d(a, infla, 0, Hhc_rot_2D, 1);

% Geometric Lower Bound
disp('==== Geometric Lower Bound ====')
[c_space3, volPolyFit] = geoLB_2d(a, infla, 0, 1e3, 1);

% "Actual" KC C-space
disp('==== "Actual" KC C-space ====')
[c_space, volMink] = mink_2d(a, infla, 0, 100, 1);

legend('Extreme vertices for Convex Lower Bound', 'Convex Lower Bound',...
    'Geometric Lower Bound','"Actual" KC C-space')