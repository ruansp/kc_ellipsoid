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

load('Hhc_3D.mat')
a = [4;2.5;2];
infla = 0.1;

b = a*(1+infla);
E1 = diag(a.^(-2));
E2 = diag(b.^(-2));

%% max angle around each semi-axis
w1 = maxAngle(a(2)/a(3), infla);
w2 = maxAngle(a(1)/a(3), infla);
w3 = maxAngle(a(1)/a(2), infla);

%% Construct c-space based on different methods
disp('==== Construct lower bounds ====')
% Convex Lower Bound
disp('-- Convex Lower Bound --')
[Z_extreme, volPoly] = cvxLB_3d(a, infla, Hhc_3D);

% Geometric Lower Bound
disp('-- Geometric Lower Bound --')
[c_space3, volPolyFit] = geoLB_3d(a, infla, 10);

%% Query points
disp('==== Query process ====')
num = [1,10,20,50,100,200,500,1e3,5e3,1e4];

for k = 1:size(num,2)
    disp(num2str([k,num(k)]))
    pt = [];
    pt = [a*infla; w1;w2;w3]*(2*rand(1,num(k))-1);
    for i = 1:size(pt,2)
        aVec = [a;pt(4:6,i)];
        bVec = [b;0;0;0];
        
        % Convex lower bound
        tp = tic;
        inPoly = inpoly(Z_extreme', pt(:,i));
        Tp(i) = toc(tp);
        
        % Geometric lower bound
        tpf = tic;
        inPolyFit(i) = ptInShrunkPoly(aVec, bVec, pt(1:3,i));
        Tpf(i) = toc(tpf);
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
ylabel('Running time (s)')
title('Containment Checking (3D case)')

