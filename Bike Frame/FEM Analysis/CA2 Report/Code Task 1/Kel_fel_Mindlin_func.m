function [K_e, K_uu, K_ww, K_wtheta, K_thetatheta_B, K_thetatheta_S, fe_ext] = Kel_fel_Mindlin_func(Ex, Ey, D_bar, D_0, G_0, p)

%% Components computable through a 4 points Gauss integration

% xi-coordinates Vector
xi_vect_4 = [-1/sqrt(3) -1/sqrt(3)  1/sqrt(3)  1/sqrt(3);
             -1/sqrt(3)  1/sqrt(3) -1/sqrt(3)  1/sqrt(3)];

% Integration weights
H_vect_4 = ones(1,4);

% Initialitations
K_uu = zeros(8);
K_thetatheta_B = zeros(8);
fe_ext_w = zeros(4,1);

% Cycle over the 4 integration points
for ii = 1:4
    xi = xi_vect_4(:, ii); 
    H = H_vect_4(ii);
    
    % Compute Determinant
    det_Fisop = detFisop_4node_func(xi, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);

    % Compute Be and Ne matrices
    [Be_u, Be_theta] = Be_Mindlin_func(xi, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);
    [~, ~, Ne_w] = N_Mindlin_func(xi, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);

    % Compute K matrices
    K_uu = K_uu + Be_u'*D_0*Be_u*det_Fisop*H;
    K_thetatheta_B = K_thetatheta_B + Be_theta'*D_bar*Be_theta*det_Fisop*H;

    % Compute Force
    fe_ext_w = fe_ext_w + Ne_w'*p*det_Fisop*H;
end
    
    
%% Components computable only through a 1 points Gauss integration (Due to shear locking)
% (Only 1 integration point for 4 nodes)

% xi-coordinates Vector
xi_1 = [0; 0];

% Integration Weight 
H_1 = 4;

% Compute Determinant
det_Fisop = detFisop_4node_func(xi_1, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);

% Compute Be and Ne matrices
[~, ~, Be_w] = Be_Mindlin_func(xi_1, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);    
[~, Ne_theta] = N_Mindlin_func(xi_1, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);

% Compute K matrices
K_ww = Be_w'*G_0*Be_w*det_Fisop*H_1;
K_wtheta = - Be_w'*G_0*Ne_theta*det_Fisop*H_1;
K_thetaw = K_wtheta';
K_thetatheta_S = Ne_theta'*G_0*Ne_theta*det_Fisop*H_1;

%% Assemblying Ke
K_e = zeros(20);

K_e(1:8, 1:8) = K_uu;
K_e(9:12, 9:12) = K_ww;
K_e(9:12, 13:20) = K_wtheta;
K_e(13:20, 9:12) = K_thetaw;
K_e(13:20, 13:20) = K_thetatheta_S + K_thetatheta_B;

% Reordering
K_e = K_e([1 2 9 13 14 3 4 10 15 16 5 6 11 17 18 7 8 12 19 20], [1 2 9 13 14 3 4 10 15 16 5 6 11 17 18 7 8 12 19 20]);

%% Assemblying fe_ext
fe_ext = zeros(20, 1);

% Inserting Forces (all the other kind of forces are 0 here)
fe_ext(9:12) = fe_ext_w;

% Reordering
fe_ext = fe_ext([1 2 9 13 14 3 4 10 15 16 5 6 11 17 18 7 8 12 19 20]');


end