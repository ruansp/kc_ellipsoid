function vol = vol_poly(poly)
%Vol_poly computes the volume of a polyhedron in R^n
%  Input:
%    poly: vertex points that defines the polyhedron in R^n
%  Output:
%    vol: Volume of the polyhedron
%
%  Author: Sipu Ruan, ruansp@jhu.edu, 2019

vol = 0;

% get simplex decomposition
enum = delaunayn(poly', {'Qt','Qbb','Qc','Qz'});

for i = 1:size(enum,1)
    simp = poly(:,enum(i,:));   
    vol = vol + vol_simplex(simp);
end

end

