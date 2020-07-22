function [A, B] = calc_derivs(x,u)
  eps = 0.0001;
  x1 = x + ones(6,1) * eps;
  x2 = x - ones(6,1) * eps;
                
  f1 = plant_dynamics(x1, u);
  f2 = plant_dynamics(x2, u);
  A = (f1 - f2) / 2 / eps;
   
  u1 = u + ones(3,1) * eps;
  u2 = u - ones(3,1) * eps;
  f3 = plant_dynamics(x, u1);
  f4 = plant_dynamics(x, u2);
  B = (f3 - f4) / 2 / eps;
end