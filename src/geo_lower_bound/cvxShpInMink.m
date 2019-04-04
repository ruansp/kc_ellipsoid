function [X_fit_poly, area] =  cvxShpInMink( aVec, bVec, N)
%MINKSUMELLIP Computes closed-form Minkowski sum and difference btw
%ellipsoids, computes convex shape inscribed inside and their areas
% Inputs:
%   aVec, bVec : vectors storing semi-axis length and rotational angle,
%                Ea is smaller than Eb,
%                ie: 2D: aVec = [a1,a2,theta_a],
%                    3D: aVec = [a1,a2,a3,w1,w2,w3];
%   N          : Number of discrete points;
%
% Output:
%   X_fit_poly : Polyhedron fit inside Minkowski difference boundary
%   area       : Area of the polyhedron
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

%% Parameters
% 2D case
if length(aVec) == 3
    dim = 2;
    phi = 0:2*pi/(N-1):2*pi;
    a(1) = aVec(1); a(2) = aVec(2); tha = aVec(3);
    b(1) = bVec(1); b(2) = bVec(2); thb = bVec(3);

    Ra = rot2(tha);
    Rb = rot2(thb);
    
% 3D case
elseif length(aVec) == 6
    dim = 3;
    a = aVec(1:3);
    b = bVec(1:3);
    
    Ra = expm(skew(aVec(4:6)));
    Rb = expm(skew(bVec(4:6)));
end

r = min(a);
T = Ra*diag(a/r)*Ra';
Binv2 = Rb*diag(1./b.^2)*Rb';

%% Shrunk space
% Find largest semi-axis
[V,D] = eig(T*Binv2*T);
[DD, DD_idx] = sort(diag(D).^(-1/2), 'descend');
VV = V(:,DD_idx);

% Extreme distance along each semi-axis
d_ext = extremeDistAxis(DD,r);

X_fit_poly = []; area = 0;
if isempty(d_ext), return; end

% Fit polyhedron
ex_vtx = diag(d_ext);
X_fit_vtx = [ex_vtx, -ex_vtx];

%% Tranform back the polyhedron, compute areas
X_fit_poly = T*VV*X_fit_vtx;

area = det(diag(a/r)) *2^dim/factorial(dim)*det(ex_vtx);

end

