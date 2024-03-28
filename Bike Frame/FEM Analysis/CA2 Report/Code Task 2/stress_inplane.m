function [sigma] = stress_inplane (Ex, Ey, ae_u, D_Mat)

% xi-coordinate of the midpoint
xi = [0; 0];

% Computing B_u matrix
Be_u = Be_Mindlin_func(xi, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);

% Computing the stresses
sigma = D_Mat*Be_u*ae_u;

end