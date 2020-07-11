function [transMatrix] = dh2transMatrix(theta, d, alpha, a)
  % Create a function that takes a set of DH parameters (each of the 4 parameters as an nx1 vector) for a link and
  % returns a 4x4 homogeneous transformation matrix
rotOldZAxis = [cos(theta) -sin(theta) 0 0;
               sin(theta) cos(theta) 0 0;
               0 0 1 0;
               0 0 0 1];
translationOldZAxis = [1 0 0 0;0 1 0 0;0 0 1 d;0 0 0 1];
translationNewXAxis = [1 0 0 a;0 1 0 0;0 0 1 0;0 0 0 1];
rotNewXAxis = [1 0 0 0;
               0 cos(alpha) -sin(alpha) 0;
               0 sin(alpha) cos(alpha) 0;
               0 0 0 1];       
transMatrix = rotOldZAxis* translationOldZAxis* translationNewXAxis* rotNewXAxis;
end