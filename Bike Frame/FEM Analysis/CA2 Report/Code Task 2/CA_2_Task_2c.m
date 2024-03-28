clc
clear
close all

%% %%%%%%%%%%%%%%% Computer Assignment 2 - Task_2 c) %%%%%%%%%%%%%%%%%%%%

% Loading Data
data_out = DataGen_CA2_Task2;

L = data_out.L;
W = data_out.W;
h = data_out.h;
D_Mat = data_out.D_Mat;
Compr_Strain = data_out.Compr_Strain;
Compr_Strenght = data_out.Compr_Strenght;

% Material Matrices computations
D_bar = h^3/12*D_Mat;
D_0 = h*D_Mat;


%% Solution of the in-plane deformation

% Creating the mesh
n_el_x = 50;
n_el_y = 5;
n_type = 2;             % 2 In-Plane solution 
[mesh, Coord, Edof] = rectMesh(0, L, 0, W, n_el_x, n_el_y, n_type);

% Element nodes Coordinates
Coord_x = Coord(:, 1);
Edof_x = Edof(:, 2:n_type:end);
Ex = Coord_x(ceil(Edof_x/n_type));

Coord_y = Coord(:, 2);
Edof_y = Edof(:, 3:n_type:end);
Ey = Coord_y(ceil(Edof_y/n_type));

% Plot mesh
figure(1)
plotpar = [1 1 0];
eldraw2(Ex,Ey,plotpar); grid on; hold on;
title('Deformation of the Specimen with prescribed displacement', 'FontSize',15);
xlabel('[mm]', 'FontSize',15); ylabel('[mm]', 'FontSize',15)

% Numbers
n_el = size(mesh, 2);
n_nodes = size(Coord, 1);
n_dof = n_nodes*n_type;

% Picking the contrained nodes
Nodes_Left = find(Coord(:, 1) == 0);
Nodes_Right = find(Coord(:, 1) == L);

% Applying Contrains 
Dof = reshape(1:n_dof, n_type, [])';
Dof_Constr_x = [Dof(Nodes_Right,1); Dof(Nodes_Left,1)];
Dof_Constr_y = [Dof(Nodes_Right,2); Dof(Nodes_Left,2)];

Dof_Constr = unique([Dof_Constr_x; Dof_Constr_y]);
n_constr = length(Dof_Constr);

Dof_Free = 1:n_dof;
Dof_Free(Dof_Constr) = [];

% Displacements
u = Compr_Strain*L;
a_uu = zeros(n_dof, 1);
a_constr = zeros(n_dof, 1);
a_constr(Dof(Nodes_Left,1)) = u;
a_constr = a_constr(Dof_Constr);

% Initializations
K_uu = spalloc(n_dof,n_dof,20*n_dof);
f_ext_uu = zeros(n_dof,1);

% Cycle over the elements
for el = 1:n_el
    [K_uu_el, ~, ~, f_el_ext] = Kirch_Quad_Routine(Ex(el,:), Ey(el, :), D_Mat, zeros(2), h, 0);
    [K_uu, f_ext_uu] = assem(Edof(el, :), K_uu, K_uu_el, f_ext_uu, f_el_ext);
end

% Solving system of equations
a_free_uu = K_uu(Dof_Free, Dof_Free)\(f_ext_uu(Dof_Free) - K_uu(Dof_Free, Dof_Constr)*a_constr);
a_uu(Dof_Free) = a_free_uu;
a_uu(Dof_Constr) = a_constr;

% Plotting Defromation
Ed = extract_dofs(Edof, a_uu);
eldisp2(Ex, Ey, Ed, [2 4 0], 1); axis tight; axis equal

figure('WindowState', 'maximized')
subplot(311)
Ed_x_avg = sum(Ed(:,1:2:end), 2)/4;
fill(Ex', Ey', Ed_x_avg); colorbar; axis equal; axis tight
title('Averaged Displacement along x   [mm]', 'FontSize',15)
xlabel('[mm]', 'FontSize',15); ylabel('[mm]', 'FontSize',15)



%% In-plane Stresses

sigma_xx = zeros(1, n_el);
N_sec = zeros(2,2, n_el);

for el = 1:n_el
    % Extraction of the element displacement
    ae_u = a_uu(Edof(el, 2:end));

    % Computing the stress
    [sigma] = stress_inplane(Ex(el, :), Ey(el, :), ae_u, D_Mat);

    % Allocating stresses
    sigma_xx(el) = sigma(1);
    N_sec(:, :, el) = h*[sigma(1) sigma(3);
                         sigma(3) sigma(2)];

end

% Plotting
figure(2)
subplot(312)
fill(Ex', Ey', sigma_xx); colorbar; axis equal; axis tight
title('Stress along x  \sigma_{xx} [MPa]', 'FontSize',15)
xlabel('[mm]', 'FontSize',15); ylabel('[mm]', 'FontSize',15)
subplot(313)
buckl = zeros(1, n_el);
buckl(:, sigma_xx < Compr_Strenght) = ones(1, size(nonzeros(sigma_xx < Compr_Strenght), 2));
fill(Ex', Ey', buckl); colorbar; axis equal; axis tight
title('Elements at Risk of failing', 'FontSize',15)
xlabel('[mm]', 'FontSize',15); ylabel('[mm]', 'FontSize',15)


%% New Boundary Conditions

% Creating the mesh
n_el_x = 50;
n_el_y = 5;
n_type = 3;             % 3 To compute the buckling 
[mesh, Coord, Edof] = rectMesh(0, L, 0, W, n_el_x, n_el_y, n_type);

% Element nodes Coordinates
Coord_x = Coord(:, 1);
Edof_x = Edof(:, 2:n_type:end);
Ex = Coord_x(ceil(Edof_x/n_type));

Coord_y = Coord(:, 2);
Edof_y = Edof(:, 3:n_type:end);
Ey = Coord_y(ceil(Edof_y/n_type));

% Numbers
n_el = size(mesh, 2);
n_nodes = size(Coord, 1);
n_dof = n_nodes*n_type;

% Picking the contrained nodes
Nodes_Left = find(Coord(:, 1) == 0);
Nodes_Right = find(Coord(:, 1) == L);
Nodes_Middle = find(Coord(:, 1) == L/2);

% Applying Contrains 
Dof = reshape(1:n_dof, n_type, [])';
Dof_Constr_w = [Dof(Nodes_Right,1); Dof(Nodes_Left,1); Dof(Nodes_Middle,1)];
Dof_Constr_omega_x = [Dof(Nodes_Right,2); Dof(Nodes_Left,2)];
Dof_Constr_omega_y = [Dof(Nodes_Right,3); Dof(Nodes_Left,3)];
Dof_Constr = unique([Dof_Constr_w; Dof_Constr_omega_x; Dof_Constr_omega_y]);
n_constr = length(Dof_Constr);


% Extraction of the free dof
Dof_Free = 1:n_dof;
Dof_Free(Dof_Constr) = [];

% Initializations
K_ww_K = spalloc(n_dof,n_dof,20*n_dof);
G_R = spalloc(n_dof,n_dof,20*n_dof);
f_ext_ww = zeros(n_dof,1);

% Cycle over the elements
for el = 1:n_el
    [~, K_ww_K_el, G_R_el] = Kirch_Quad_Routine(Ex(el,:), Ey(el, :), D_Mat, N_sec(:,:, el), h, 0);
    K_ww_K = assem(Edof(el, :), K_ww_K, K_ww_K_el);
    G_R = assem(Edof(el, :), G_R, G_R_el);
end



%% Eigenvalues problem

% Extracting Free part of the matrices
K_ww_K_FF = K_ww_K(Dof_Free, Dof_Free);
G_R_FF = G_R(Dof_Free, Dof_Free);

% Solving Eigenvalues problems
n_lambda = 8;
[Eig_Vect, Eig_Val_Mat] = eigs(K_ww_K_FF, -G_R_FF, n_lambda, 'smallestabs');
lambda = diag(Eig_Val_Mat);

% Displaying Safety Factor
fprintf(['The Safety Factor against buckling is: ', num2str(lambda(1))]);

figure('WindowState', 'maximized')
sgtitle('Mode Shapes', 'FontSize', 18, 'FontWeight', 'bold')
for ii = 1:n_lambda
    % Find the modeshape
    Mode_Shape = zeros(n_dof, 1);
    Mode_Shape(Dof_Free, :) = Eig_Vect(:, ii);

    % Extraction to plot
    E_Mode_Shape = extract_dofs(Edof, Mode_Shape);
    E_Mode_Shape_avg = sum(E_Mode_Shape(:,1:3:end), 2)/4;

    % Plotting
    subplot(ceil(n_lambda/2), 2, ii)
    fill(Ex', Ey', E_Mode_Shape_avg); colorbar; axis equal; axis tight
    title(['Mode Shape Number: ', num2str(ii), '      \lambda = ', num2str(lambda(ii))], 'FontSize',15); 
    xlabel('[mm]', 'FontSize',15); ylabel('[mm]', 'FontSize',15)
end
