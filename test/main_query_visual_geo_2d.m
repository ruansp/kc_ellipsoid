% Visualization of Geometric Lower Bound (2D) and its containment checking
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
alpha = a(1)/a(2);
infla = .1;
b = a*(1+infla);
A = [a(1)^(-2) 0;0 a(2)^(-2)];
VA = [a(1) 0; 0 a(2)];
B = [b(1)^(-2) 0;0 b(2)^(-2)];

% Max angle of rotation
phi = maxAngle(alpha,infla);

%% Geometric Lower Bound
disp('== Computing Geometric Lower Bound ==')
[c_space, ~] = geoLB_2d(a, infla, 0, 1e3, 0);

figure(1); hold on; axis equal; grid on;
plot3(c_space(1,:),c_space(2,:),c_space(3,:),'-.')

%% Query points in c-space
disp('== Containment Checking ==')
num = 1e3;
pt = [(2*rand(1,num)-1)*0.6;...
    (2*rand(1,num)-1)*0.4;...
    (2*rand(1,num)-1)*phi];
lw = 1;

figure(2); hold on; axis equal; grid on;
% Larger ellipse
bVec = [b;0];
Ellipse_plot(diag(b)^(-2), [0;0], 'black');
for i = 1:size(pt,2)
    aVec = [a;pt(3,i)];
    
    in(i) = ptInShrunkPoly(aVec, bVec, pt(1:2,i));
    
    if in(i)
        figure(1);
        plot3(pt(1,i),pt(2,i),pt(3,i),'g+','LineWidth',lw)
        
        figure(2);
        R = rot2(pt(3,i));
        A = R*diag(a)^(-2)*R';
        Ellipse_plot(A,[pt(1,i); pt(2,i)], 'green');
    else
        figure(1);
        plot3(pt(1,i),pt(2,i),pt(3,i),'r.','LineWidth',lw)
    end
end

% Check algebraic condition
pnt_valid = [];
ang = 0:pi/20:2*pi;
u = [cos(ang); sin(ang)];
for i = 1:size(pt,2)
    if in(i)
        inside = 1;
        for j = 1:size(ang,2)
            % condition that smaller ellipse moves inside larger one without
            % collision.
            R = expm(skew(pt(3,i)));
            t = pt(1:2,i);
            ineqF = (R*VA*u(:,j) + t')' * B * (R*VA*u(:,j) + t');
            if ineqF - 1 > 1e-8
                inside = 0;
                break;
            end
        end
        if inside
            pnt_valid = [pnt_valid; pt(:,i)'];
        end
    end
end
disp('Configs inside Convex Lower Bound: ')
pt(:,in==1)'
disp('Valid configs from Algebraic Condition of Containment: ')
pnt_valid
