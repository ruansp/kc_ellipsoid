function R = rot2(th)
%ROT2 Rotation matrix in 2D, i.e. SO(2)

R = [cos(th), -sin(th); 
     sin(th),  cos(th)];
end

