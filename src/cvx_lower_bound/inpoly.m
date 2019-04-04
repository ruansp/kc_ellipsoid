function in = inpoly( poly, pnt )
%INPOLY Checks whether a point in R^n is in the interior of the polyhedron
%with (n+1) vertices.
% Inputs:
%   poly:    the polyhedron vertices in R^n
%   pnt:     point to be checked
% Output:
%   in:      result
% Example:
%   poly = [1, 1, -1, -1, 1; 1, -1, -1, 1, -3; 1, -1, 1, -1, 5];
%   pnt = rand(3,1);
%   in = inpoly(poly, pnt);
% 
%   K = convhull(poly(1,:), poly(2,:), poly(3,:));
%   patch('Faces',K, 'Vertices',[poly'], 'FaceColor', 'y', 'FaceAlpha', 0.3);
%   hold on;
%   if in
%       plot3(pnt(1,:), pnt(2,:), pnt(3,:), 'g*');
%   else
%       plot3(pnt(1,:), pnt(2,:), pnt(3,:), 'ro');
%   end
% 
% Author: Sipu Ruan, ruansp@jhu.edu, August 2019

% initialization
if size(poly,1) ~= size(pnt, 1)
    error('Dimensions of polyhedron vertices and point do not match! ');
end

in = zeros(1,size(pnt,2));

% get simplex decomposition
enum = delaunayn(poly');

for i = 1:size(enum,1)
    simp = poly(:,enum(i,:));
    in = in + insimplex(simp, pnt);
end

in(in>0) = 1;