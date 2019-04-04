function [c_space, volMink] = mink_3d(a, infla, N)
% mink_3d constructs KC space as a union of Minkowski difference at
% different angles
%
% Inputs:
%  a     : semi-axis length of the smaller ellipse
%  infla : inflation factor ( >0 )
%  th_b  : rotational angle of the larger ellipse
%  N     : number of sampled points in a sphere
%  isplot: whether to plot the c-space
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

%% Parameters
b = a*(1+infla);

% Initialize ET parameters
global ellOptions;
ellOptions.verbose = 0;
ellOptions.rel_tol = 1e-6;
ellOptions.abs_tol = 1e-7;
ellOptions.plot3d_grid = 150;

VA = diag(a.^2);
VB = diag(b.^2);

%% Exponential coordinates for rotation
% max angle around each semi-axis
w1 = maxAngle(a(2)/a(3), infla);
w2 = maxAngle(a(1)/a(3), infla);
w3 = maxAngle(a(1)/a(2), infla);

alpha = -w1:2*w1/(N-1):w1; 
beta = -w2:2*w2/(N-1):w2; 
gamma = -w3:2*w3/(N-1):w3;
[aa,bb,cc] = ndgrid(alpha, beta, gamma);
ww = [aa(:),bb(:),cc(:)]';

%% For each angle, compute mink diff
c_space = [];
intgrand = zeros(1,size(ww,2));
for i = 1:size(ww,2)
    % Closed-form Mink diff
    aVec = [a;ww(:,i)];
    bVec = [b;zeros(3,1)];
    Ra = expm(skew(aVec(4:6)));
    Rb = expm(skew(bVec(4:6)));
    
    % Ellipsoid using ET
    A = Ra * VA * Ra';
    B = Rb * VB * Rb';
    
    A(1,2) = A(2,1); A(1,3) = A(3,1); A(2,3) = A(3,2);
    B(1,2) = B(2,1); B(1,3) = B(3,1); B(2,3) = B(3,2);
    
    Ea = ellipsoid(A);
    Eb = ellipsoid(B);
    
    % Compute Minkowski difference
    [~, bd{i}] = minkdiff(Eb, Ea);

    % Volume and store boundary points
    if size(bd{i},2)<=4
        i = i+1;
        continue;
    end
    
    [~,V_ang(i)] = convhull(bd{i}');
    detJac(i) = 2*(1-cos(norm(ww(:,i))))/norm(ww(:,i))^2;
    intgrand(i) = V_ang(i) * detJac(i);
    
    c_space = [c_space [bd{i}; ww(:,i)*ones(1,size(bd{i},2))]];
end

intgrand = reshape(intgrand, [N,N,N]);
volMink = trapz(alpha,trapz(beta,trapz(gamma,intgrand,3),2),1);