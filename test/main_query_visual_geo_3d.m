% Geometric Lower Bound (3D) containment checking
%
% Dependency:
%  "Ellipsoidal Toolbox" for plotting ellipsoids
%  "Robotics Toolbox" by Peter Corkes for rotations
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

close all; clear; clc;

%% Parameters
addpath ../include/
addpath ../mat/
addpath ../src/cvx_lower_bound/
addpath ../src/geo_lower_bound/

a = [4;2.5;2];
infla = .1;
b = a*(1+infla);
bVec = [b;zeros(3,1)];

VA = diag(a.^2);
B = diag(b.^2);

% max angle around each semi-axis
w1 = maxAngle(a(2)/a(3), infla);
w2 = maxAngle(a(1)/a(3), infla);
w3 = maxAngle(a(1)/a(2), infla);

%% Query points in C-space
disp('== Containment Checking for Geometric Lower Bound ==')
num = 1e3;
pt = [a*infla; w1;w2;w3]*(2*rand(1,num)-1);

figure(1); hold on; axis equal; grid on;
% Larger ellipse
Eb = ellipsoid([0;0;0], B);
plot(Eb, 'b')

% Ellipse_plot(diag(b)^(-2), [0;0;0], 'b', 'none', 0.3);
for i = 1:size(pt,2)
    aVec = [a;pt(4:6,i)];
    
    in(i) = ptInShrunkPoly(aVec, bVec, pt(1:3,i));
    
    R = expm(skew(pt(4:6,i)));
    A = R*VA*R';
    A(1,2) = A(2,1); A(1,3) = A(3,1); A(2,3) = A(3,2);
    
    if in(i)
        Ea = ellipsoid([pt(4,i); pt(5,i); pt(6,i)], A);
        plot(Ea, 'g')
    end
end

% Check algebraic condition
pnt_valid = configValidation3D(pt(:,in==1)', a ,infla);
disp('Configs inside Geometric Lower Bound: ')
pt(:,in==1)'
disp('Valid configs from Algebraic Condition of Containment: ')
pnt_valid
