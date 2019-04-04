% Visualization of Convex Lower Bound (2D) and its containment checking
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

close all; clear; clc;

%% Parameters
addpath ../include/
addpath ../mat/
addpath ../src/cvx_lower_bound/
addpath ../src/geo_lower_bound/

load('Hhc_rot_2D.mat')
a = [5;3.5];
ratio = a(1)/a(2);
infla = 0.08;
b = a * (1+infla);
A = [a(1)^(-2) 0;0 a(2)^(-2)];
VA = [a(1) 0; 0 a(2)];
B = [b(1)^(-2) 0;0 b(2)^(-2)];

%% Convex Lower Bound
disp('== Computing Convex Lower Bound ==')
[Z_extreme, volCvx] = cvxLB_2d(a, infla, 0, Hhc_rot_2D, 0);
x = Z_extreme(:,1); y = Z_extreme(:,2); z = Z_extreme(:,3);

%% Test for point in polyhedron
disp('== Containment Checking ==')
numTest = 1e3;
testPnt = 2*[(rand(numTest,1)-0.5)*Z_extreme(5,1),...
    (rand(numTest,1)-0.5)*Z_extreme(7,2),...
    (rand(numTest,1)-0.5)*Z_extreme(9,3)];
in = inpoly(Z_extreme', testPnt');

% Check algebraic condition
pnt_valid = [];
ang = 0:pi/20:2*pi;
u = [cos(ang); sin(ang)];
for i = 1:size(testPnt,1)
    if in(i)
        inside = 1;
        for j = 1:size(ang,2)
            % condition that smaller ellipse moves inside larger one without
            % collision.
            R = expm(skew(testPnt(i,3)));
            t = testPnt(i,1:2);
            ineqF = (R*VA*u(:,j) + t')' * B * (R*VA*u(:,j) + t');
            if ineqF - 1 > 1e-8
                inside = 0;
                break;
            end
        end
        if inside
            pnt_valid = [pnt_valid; testPnt(i,:)];
        end
    end
end
disp('Configs inside Convex Lower Bound: ')
testPnt(in==1,:)
disp('Valid configs from Algebraic Condition of Containment: ')
pnt_valid

%% Plots
disp('== Plotting results ==')
% C-space
figure(1); hold on; axis equal; grid on;
plot3(x,y,z,'b*', 'LineWidth', 2);
[K, volPoly] = convhull(x,y,z);
P = patch('Faces',K, 'Vertices',[x y z], 'FaceColor', 'y', 'FaceAlpha', 0.3);

% Euclidean space
figure(2); hold on; axis equal; grid on;
% Larger ellipse
Ellipse_plot(B, [0;0], 'black');

% Smaller ellipse
for i = 1:size(z,1)
    R_A = rot2(z(i));
    A2 = R_A*A*R_A';
    Ellipse_plot(A2, [x(i); y(i)], 'blue');
end

% Test point-in-poly
figure(3); hold on; axis equal;
P = patch('Faces',K, 'Vertices',[x y z], 'FaceColor', 'y', 'FaceAlpha', 0.3);

figure(4); hold on; axis equal;
Ellipse_plot(B, [0;0], 'black');
num = 1;
for i = 1:size(testPnt,1)
    if in(i)
        figure(3);
        plot3(testPnt(i,1),testPnt(i,2),testPnt(i,3), 'g+', 'LineWidth', 1);

        figure(4); axis equal; hold on;
        R_A = eye(2)+skew(testPnt(i,3));
        A2 = R_A*A*R_A';
        Ellipse_plot(A2, [testPnt(i,1); testPnt(i,2)], 'green');
        num = num+1;
    else
        figure(3);
        plot3(testPnt(i,1),testPnt(i,2),testPnt(i,3), 'r.', 'LineWidth', 1);
    end
end

figure; hold on; 
plot3(pnt_valid(:,1),pnt_valid(:,2),pnt_valid(:,3),'g+')
plot3(pnt_valid2(:,1),pnt_valid2(:,2),pnt_valid2(:,3),'mo')