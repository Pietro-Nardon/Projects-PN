clc
clear 
close all

%% %%%%%%%%%%%%%%% Computer Assignment 2 - Task_1 b) %%%%%%%%%%%%%%%%%%%%%%
% Geometry
Ex_Test = [2 11 12 3]*1e-3;        % [m]     
Ey_Test = [4 5 21 22]*1e-3;        % [m]
h_test = 50e-3;                    % [m]


% Material
E_Test = 80e9;                     % [Pa]
nu_Test = 0.2;
G_modulus_Test = E_Test/(2*(1+nu_Test));
G_Test = G_modulus_Test*eye(2,2);
G_0_Test = G_Test*h_test*5/6;
D_Test = hooke(1,E_Test,nu_Test);  % Note: Thin plate = plane stress 2D assumption
D_bar_Test = h_test^3/12*D_Test;
D0_Test = D_Test*h_test;

% Displacements
ae_Test = (1:20)'*1e-6;
z_Test = 5e-3;

% Out-of-plane load qz
q_Test = 0.5*1e3;                  

% Plot Mesh
figure
plotpar = [1 1 0];
eldraw2(Ex_Test,Ey_Test,plotpar); grid on
title('Undeformed Plate Test', 'FontSize',15); 
xlabel('[m]', 'FontSize',15); ylabel('[m]', 'FontSize',15)

% Stiffness Matrix and Load Vector
[K_e_Test, K_uu_Test, K_ww_Test, K_wtheta_Test, K_thetatheta_B_Test, K_thetatheta_S_Test, fe_ext_Test] = Kel_fel_Mindlin_func(Ex_Test, Ey_Test, D_bar_Test, D0_Test, G_0_Test, q_Test);

% Stresses
[sigma_el_Test, tau_el_Test] = Stresses_Mindlin_func(Ex_Test, Ey_Test, ae_Test, z_Test, D_Test, G_Test);
