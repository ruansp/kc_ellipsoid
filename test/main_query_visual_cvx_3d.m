% Visualization of Convex Lower Bound (3D) and its containment checking
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
load('Hhc_3D.mat');

a = [4;2.5;2];
infla = 0.1;
b = a * (1+infla);
VA = diag(a.^2);
B = diag(b.^2);

%% Convex Lower Bound
disp('== Computing Convex Lower Bound ==')
[Z_extreme, ~] = cvxLB_3d(a, infla, Hhc_3D);

th1 = Z_extreme(:,1); th2 = Z_extreme(:,2); th3 = Z_extreme(:,3);
x = Z_extreme(:,4); y = Z_extreme(:,5); z = Z_extreme(:,6);

%% Test points inside the polyhedron
disp('== Containment Checking ==')
numTest = 1e3;
testPnt = [(max(th1)-min(th1))*rand(numTest,1)+min(th1)...
           (max(th2)-min(th2))*rand(numTest,1)+min(th2)...
           (max(th3)-min(th3))*rand(numTest,1)+min(th3)...
           (max(x)-min(x))*rand(numTest,1)+min(x)...
           (max(y)-min(y))*rand(numTest,1)+min(y)...
           (max(z)-min(z))*rand(numTest,1)+min(z)];
in = inpoly(Z_extreme', testPnt');

% Check algebraic condition
pnt_valid = configValidation3D(testPnt(in==1,:), a ,infla);
disp('Configs inside Convex Lower Bound: ')
testPnt(in==1,:)
disp('Valid configs from Algebraic Condition of Containment: ')
pnt_valid

%% Plots
%% C-space
disp('== Plotting results ==')
% Translation
figure(1); 
subplot(1,2,1); hold on;
title('Extreme configurations in C-space: Translation part')

plot3(x,y,z,'b*', 'LineWidth', 2);
K = convhull(x,y,z);
P = patch('Faces',K, 'Vertices',[x y z], 'FaceColor', 'y', 'FaceAlpha', 0.3);

% Rotation
subplot(1,2,2); hold on;
title('Extreme configurations in C-space: Rotation part')

plot3(th1,th2,th3,'b*', 'LineWidth', 2);
K = convhull(th1,th2,th3);
P = patch('Faces',K, 'Vertices',[th1 th2 th3], 'FaceColor', 'y', 'FaceAlpha', 0.3);

%% Euclidean space
figure(2); 
hold on; axis equal;
% Larger ellipsoid
Eb = ellipsoid([0;0;0], B);
plot(Eb, 'b')

% Smaller ellipsoid
for i = 1:size(z,1)
    R_A = expm(skew([th1(i);th2(i);th3(i)]));
    A = R_A*VA*R_A';
    A(1,2) = A(2,1); A(1,3) = A(3,1); A(2,3) = A(3,2);
    
    Ea = ellipsoid([x(i); y(i); z(i)], A);
    plot(Ea, 'g')
end
title('Exteme configurations in Euclidean space')

%% Test pnt-in-poly
figure(3); hold on; axis equal; grid on;
% Larger ellipsoid
Eb = ellipsoid([0;0;0], B);
plot(Eb, 'b')

% Smaller ellipsoid
for i = 1:size(testPnt,1)
    R_A = expm(skew([testPnt(i,1); testPnt(i,2); testPnt(i,3)]));
    A = R_A*VA*R_A';
    A(1,2) = A(2,1); A(1,3) = A(3,1); A(2,3) = A(3,2);
    
    if in(i)
        Ea = ellipsoid([testPnt(i,4); testPnt(i,5); testPnt(i,6)], A);
        plot(Ea, 'g')
    end
end
title('Result for containment test in Euclidean space')
