function [K_el, f_el_ext] = Kirch_K_f_4noded(Ex, Ey, p, D_bar)

% xi-coordinates Vector
xi_vect = [-1/sqrt(3) -1/sqrt(3)  1/sqrt(3)  1/sqrt(3);
           -1/sqrt(3)  1/sqrt(3) -1/sqrt(3)  1/sqrt(3)];

% Integration weights
H_vect = ones(1,4);

% Initializing
K_el = zeros(12,12);
f_el_ext = zeros(12,1);

% Loop over integration points
for ii = 1:4 
    xi = xi_vect(:, ii); 
    H = H_vect(ii);
    
    % Compute Determinant
    det_Fisop = detFisop_4node_func(xi, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);

    % Compute Shape Functions
    N = N_kirchoff_func(xi, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);

    % Compute element external force
    f_el_ext = f_el_ext + N'*p*det_Fisop*H;

    % Compute Bstar
    B_Star = Bast_kirchoff_func(xi, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);

    % Compute Stiffnes matrix
    K_el = K_el + B_Star'*D_bar*B_Star*det_Fisop*H;
end
end