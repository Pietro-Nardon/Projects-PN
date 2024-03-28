function [K_e, f_e] = CST_stiff_load(node_coord, D_Mat, F_g, h)
% This function computes the Element Stiffness Matrix and the Element Load Vector loaded with gravity
%
% ------------------------------------------------------
% Input:
%           node_coord = [2x3] (1 row = x; 2 row = y)
%                    E = Elastic Module of the material
%                   nu = Poisson's Ratio of the material
%                  F_g = Gravitational Force 
%                    h = Thickness of the Element
% ------------------------------------------------------
% Output:
%                  K_e = Stiffness matrix of the CST element
%                  f_e = Element Volume Load Vector


% Number of elements
n_el = 1;       % We work with CST, the nodes are just 3, so we'll have 1 triangle so 1 element


% Plot of the element
% figure
% line(node_coord(1,:), node_coord(2, :),'Color', 'red', 'LineWidth', 2); hold on
% line(node_coord(1,[3 1]), node_coord(2,[3 1]),'Color', 'red', 'LineWidth', 2); grid on;
% title('CST Element'); xlabel('[m]'); ylabel('[m]');


%% Creation of Ae Be

A_e = Ae_cst_func(node_coord(:, 1), node_coord(:, 2),node_coord(:, 3));
B_e = Be_cst_func(node_coord(:, 1), node_coord(:, 2),node_coord(:, 3));

%% Compute Ke and f

% Element stiffness matrix
K_e = B_e'*D_Mat*B_e*h*A_e;

% Body Force (due to only gravity)
F_Body = zeros(2,n_el);
F_Body(2,:) = F_g;             % Because it acts only along y

% Load vector 
f_e = [F_Body; F_Body; F_Body]*A_e/3;

end