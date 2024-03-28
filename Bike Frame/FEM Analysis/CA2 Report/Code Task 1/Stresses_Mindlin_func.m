function [sigma_el, tau_el] = Stresses_Mindlin_func(Ex, Ey, ae, z, D_Mat, G_Mat)
% This function computes the stresses in a mindlin plate element
%assuming isotropic linear elasticity with plane stress 
% Output:    sigma_el = [11 12; 21 22] tau_el = [13 23]

% xi-coordinates Vector
xi_vect = [-1/sqrt(3) -1/sqrt(3)  1/sqrt(3)  1/sqrt(3);
           -1/sqrt(3)  1/sqrt(3) -1/sqrt(3)  1/sqrt(3)];

% Displacements Extracions
ae_u = ae([1 2 6 7 11 12 16 17]);
ae_w = ae([3 8 13 18]);
ae_theta = ae([4 5 9 10 14 15 19 20]);

% Initializations
sigma_el = zeros(3, 4);
tau_el = zeros(2, 4);

% Cycle over the 4 quadrature points
for ii = 1:4
    xi = xi_vect(:, ii); 
    
    % Compute Be and Ne matrices
    [Be_u, Be_theta, Be_w] = Be_Mindlin_func(xi, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);
    [~, Ne_theta] = N_Mindlin_func(xi, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);

    % Compute stresses
    sigma_el(:, ii) = D_Mat*(Be_u*ae_u - z*Be_theta*ae_theta);
    tau_el(:, ii) = G_Mat*(Be_w*ae_w - Ne_theta*ae_theta);
end
end