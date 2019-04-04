function in = insimplex( simp, pnt )
%INSIMPLEX Checks whether a point in R^n is in the interior of the simplex
%with (n+1) vertices.
% Inputs:
%   simp:    the simplex vertices in R^n
%   pnt:     point to be checked
% Output:
%   in:      result
% Example:
%   simp = [1, 1, -1, -1; 1, -1, -1, 1; 1, -1, 1, -1];
%   pnt = rand(3,1);
%   in = insimplex(simp, pnt);
%
%   K = convhull(simp(1,:), simp(2,:), simp(3,:));
%   patch('Faces',K, 'Vertices',[simp'], 'FaceColor', 'y', 'FaceAlpha', 0.3);
%   hold on;
%   plot3(pnt(1,:), pnt(2,:), pnt(3,:), '*')
% 
% Author: Sipu Ruan, ruansp@jhu.edu, August 2019

% initialization
if size(simp,1) ~= (size(simp,2)-1)
    error('The given shape is not a simplex in R^n! ');
end

if size(simp,1) ~= size(pnt, 1)
    error('Dimensions of simplex vertices and point do not match! ');
end

S = [simp; ones(1,size(simp,2))];
P = [pnt; ones(1,size(pnt,2))];

% Calculate parameters lambda
lambda = S\P;

% Check if all lambdas are in [0,1]
in = all(lambda>=0, 1) & all(lambda<=1, 1);

end

