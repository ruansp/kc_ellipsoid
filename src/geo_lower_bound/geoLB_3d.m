function [c_space, volPolyFit] = geoLB_3d(a, infla, N)
% geoLB_3d constructs KC space as a union of closed-form fitted polyhedron
% at different orientations (Geometric Lower Bound)
%
% Inputs:
%  a     : semi-axis length of the smaller ellipse
%  infla : inflation factor ( >0 )
%  N     : number of sampled points on the surface
% 
% Outputs:
%  c_space   : vertices of configurations in C-space
%  volPolyFit: Volume of the fitted polyhedron
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

%% Parameters
b = a*(1+infla);

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
    [X_fit_poly, V_poly(i)] = cvxShpInMink(aVec, bVec, 50);
    
    % Volume and store boundary points
    if isempty(X_fit_poly)
        i = i+1;
        continue;
    end
    
    detJac(i) = 2*(1-cos(norm(ww(:,i))))/norm(ww(:,i))^2;
    intgrand(i) = V_poly(i) * detJac(i);
    
    c_space = [c_space [X_fit_poly; ww(:,i)*ones(1,size(X_fit_poly,2))]];
end

intgrand = reshape(intgrand, [N,N,N]);
volPolyFit = trapz(alpha,trapz(beta,trapz(gamma,intgrand,3),2),1);