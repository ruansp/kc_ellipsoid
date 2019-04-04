function d_ext = extremeDistAxis(a,r)
% extrameDistAxis computes the extreme distance a sphere can travel at the
% semi-axis of an ellipsoid while fully contained.
%  Inputs:
%    a: semi-axis lengths
%    r: radius of the sphere
%  Outputs:
%    d: extreme distance
%
%  Author: Sipu Ruan, ruansp@jhu.edu, 2019

% tic;
d_ext = [];
dim = length(a);
for i = 1:dim
    if r > a(i)
%         warning('Sphere radius should not be larger than any semi-axis length')
        return;
    end
    
    for j = 1:dim
        if r >= a(j)^2/a(i)
            dd(i,j) = sqrt( (a(j)^2-a(i)^2)*(r^2-a(j)^2)/a(j)^2 );
        else
            dd(i,j) = a(i) - r;
        end
    end
end
d_ext = min(dd,[],2);