% Comparisons of containment checking time between different methods in 3D
% 
% Variables:
%   Random configurations
%   Number of query points
% 
% C-space to be compared:
%  (1) Convex Lower Bound
%  (2) Geometric Lower Bound
%
% Plots:
%  (1) Running time with respect to the number of queries
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

clc; clear; close all;

%% Parameters
addpath ../include/
addpath ../mat/
addpath ../src/cvx_lower_bound/
addpath ../src/geo_lower_bound/

load('Hhc_rot_2D.mat')
a = [5;3.5];
infla = 0.1;
alpha = a(1)/a(2);

b = a*(1+infla);
E1 = diag(a.^(-2));
E2 = diag(b.^(-2));

%% Max angle of rotation
w = maxAngle(alpha,infla);

%% Construct c-space based on different methods
disp('==== Construct lower bounds ====')
% Convex Lower Bound
disp('-- Convex Lower Bound --')
[Z_extreme, volPoly] = cvxLB_2d(a, infla, 0, Hhc_rot_2D, 0);

% Geometric Lower Bound
disp('-- Geometric Lower Bound --')
[c_space3, volPolyFit] = geoLB_2d(a, infla, 0, 1e3, 0);

%% Query points
disp('==== Query process ====')
num = [1,10,20,50,100,200,500,1e3,5e3,1e4];

for k = 1:size(num,2)
    disp(num2str([k,num(k)]))
    pt = [];
    pt = [a*infla; w]*(2*rand(1,num(k))-1);
    for i = 1:size(pt,2)
        aVec = [a;pt(3,i)];
        bVec = [b;0];
        
        % Geometric lower bound
        tpf = tic;
        inPolyFit(i) = ptInShrunkPoly(aVec, bVec, pt(1:2,i));
        Tpf(i) = toc(tpf);
        
        % Convex lower bound
        tp = tic;
        inPoly = inpoly(Z_extreme', pt(:,i));
        Tp(i) = toc(tp);
    end
    
    mTpf(k) = sum(Tpf);
    mTp(k) = sum(Tp);
end


%% Plots
figure; hold on; grid on;
lw = 1.25;
plot(num, mTp, 'b-.', 'LineWidth', lw);
plot(num, mTpf, 'r--', 'LineWidth', lw);

legend('Convex Lower Bound', 'Geometric Lower Bound')
xlabel('Number of sample points')
ylabel('Running time (s)')
title('Containment Checking (2D)')

