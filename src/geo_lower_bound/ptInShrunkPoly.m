function in = ptInShrunkPoly(aVec, bVec, pt)
% ptInShrunkPoly checks whether a given point is inside the polyhderon in
% shrunk space
%
% Inputs:
%  aVec, bVec: a structured the vector for the shape of two ellipsoids,
%              i.e. in 2D, aVec = [a1,a2,theta];
%                   in 3D, aVec = [a1,a2,a3,w1,w2,w3];
%  pt        : queried point that defines the configuration
% 
% Outputs:
%  in        : results of the containment checking
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2019

%% Parameters
% 2D case
if length(aVec) == 3
    a(1) = aVec(1); a(2) = aVec(2); tha = aVec(3);
    b(1) = bVec(1); b(2) = bVec(2); thb = bVec(3);
    Ra = rot2(tha);
    Rb = rot2(thb);
    
    r = min(a);
    T = Ra*diag(a/r)*Ra';
    B = Rb*diag(b)*Rb';
    Binv2 = Rb*diag(b.^(-2))*Rb';

    % 3D case
elseif length(aVec) == 6
    a = aVec(1:3);
    b = bVec(1:3);
    
    Ra = expm(skew(aVec(4:6)));
    Rb = expm(skew(bVec(4:6)));
    
    r = min(a);
    T = Ra*diag(a/r)*Ra';
    B = Rb*diag(b)*Rb';
    Binv2 = Rb*diag(1./b.^2)*Rb';
end

%% Test for Point-in-Ellipsoid
[V,D] = eig(T*Binv2*T);
[DD, DD_idx] = sort(diag(D).^(-1/2), 'descend');
VV = V(:,DD_idx);

% Transform point to the shrunk ellipse coord
pt_t = VV'*(T\pt);

in = false;
d_ext = extremeDistAxis(DD,r);
if isempty(d_ext), return; end

ex_vtx = diag(d_ext);

poly = [ex_vtx, -ex_vtx];
in = inpoly(poly, pt_t);
