function vol = vol_simplex(P)
%Vol_simplex computes the volume of an n-simplex in R^n
%  Input:
%    P: n x (n+1) matrix storing n+1 vertex points in R^n
%  Output:
%    vol: Volume of the simplex
%
%  Author: Sipu Ruan, ruansp@jhu.edu, 2019

if size(P,1)+1 ~= size(P,2)
    error('Not an n-simplex in R^n, should input an n x (n+1) matrix');
end

n = size(P,1);
diffMat = P(:,2:n+1)-P(:,1);
vol = abs( 1/factorial(n) * det(diffMat) );

end

