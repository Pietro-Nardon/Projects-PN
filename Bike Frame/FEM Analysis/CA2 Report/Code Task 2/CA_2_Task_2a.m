clc
clear
close all

%% %%%%%%%%%%%%%%% Computer Assignment 2 - Task_2 a) %%%%%%%%%%%%%%%%%%%%%%

% Geometry
Ex_Test = [2 11 12 3]*1e-3;        % [m]     
Ey_Test = [4 5 21 22]*1e-3;        % [m]
h_test = 50e-3;                    % [m]

% Material
E_Test = 80e9;                     % [Pa]
nu_Test = 0.2;
D_Test = hooke(1,E_Test,nu_Test);  % Note: Thin plate = plane stress 2D assumption

% Stresses
sigma2D = [1 2; 2 3]*1e6;
N_sec = h_test*sigma2D;

% Matrices
[K_uu_el, K_ww_K_el, G_el] = Kirch_Quad_Routine(Ex_Test, Ey_Test, D_Test, N_sec, h_test, -1);

% Display
disp(K_ww_K_el); disp(G_el); disp(K_uu_el)

