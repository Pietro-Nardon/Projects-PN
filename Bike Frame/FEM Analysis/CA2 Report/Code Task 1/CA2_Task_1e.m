clc
clear
close all

%% %%%%%%%%%%%%%%% Computer Assignment 2 - Task_1 e) %%%%%%%%%%%%%%%%%%%%%%


%% Loading Datas
data = DataGen_CA2;

L = data.L;
W = data.W;
h_min = data.h_min;
h_max = data.h_max;
E_Mat = data.E_Mat;
nu = data.nu;
p = data.p;
h = 10;             %[mm]

% Material Matrix
ptype = 1;              % Plane Stress
D_Mat = hooke(ptype, E_Mat, nu);
D_bar = h^3/12*D_Mat;
D_0 = h*D_Mat;
G_Mat = E_Mat/(2*(1+nu));
G_0 = G_Mat*5/6*h;

% Creating the mesh
n_el_x = 20;
n_el_y = 30;
n_type = 5;             % 3 Kirch       % 5 Mindlin
[mesh, Coord, Edof] = rectMesh(0, W, 0, L, n_el_x, n_el_y, n_type);

% Element nodes Coordinates
Coord_x = Coord(:,1);
Edof_x = Edof(:, 2:n_type:end);
Ex = Coord_x(ceil(Edof_x/n_type));

Coord_y = Coord(:, 2);
Edof_y = Edof(:, 3:n_type:end);
Ey = Coord_y(ceil(Edof_y/n_type));

% Plot mesh
figure(1)
plotpar = [1 1 0];
eldraw2(Ex,Ey,plotpar); grid on
title('Undeformed Plate with Mindlin Elements', 'FontSize',15); axis tight
xlabel('[mm]', 'FontSize',15); ylabel('[mm]', 'FontSize',15)

% Numbers
n_el = size(mesh, 2);
n_nodes = size(Coord, 1);
n_dof = n_nodes*n_type;

% Picking dofs of the edge elements
Nodes_Bottom = find(Coord(:, 2) == 0);
Nodes_Top = find(Coord(:, 2) == L);
Nodes_Left = find(Coord(:, 1) == 0);
Nodes_Right = find(Coord(:, 1) == W);

% Applying Contrains 
Dof = reshape(1:n_dof, n_type, [])';
Dof_Constr_x = [Dof(Nodes_Bottom, 1); Dof(Nodes_Right,1); Dof(Nodes_Left,1); Dof(Nodes_Top,1)];
Dof_Constr_y = [Dof(Nodes_Bottom, 2); Dof(Nodes_Right,2); Dof(Nodes_Left,2); Dof(Nodes_Top,2)];
Dof_Constr_w = [Dof(Nodes_Bottom, 3); Dof(Nodes_Right,3); Dof(Nodes_Left,3); Dof(Nodes_Top,3)];
Dof_Constr_Rot = [Dof(Nodes_Bottom, 4); Dof(Nodes_Right,5); Dof(Nodes_Left,5); Dof(Nodes_Top,4)];
Dof_Constr = unique([Dof_Constr_x; Dof_Constr_y; Dof_Constr_w; Dof_Constr_Rot]);
n_constr = length(Dof_Constr);

Dof_Free = 1:n_dof;
Dof_Free(Dof_Constr) = [];



%% Computation of the displacements

% Initializations
K = spalloc(n_dof,n_dof,20*n_dof);
f_ext = zeros(n_dof,1);
a = zeros(n_dof, 1);
a_constr = zeros(n_constr, 1);


% Cycle looping over all the elements to compute K and f
for el = 1:n_el
    [K_el, ~, ~, ~, ~, ~, f_el_ext] = Kel_fel_Mindlin_func(Ex(el,:), Ey(el, :), D_bar, D_0, G_0, p);
    [K, f_ext] = assem(Edof(el, :), K, K_el, f_ext, f_el_ext);
end

% Solving system of equations
a_free = K(Dof_Free, Dof_Free)\(f_ext(Dof_Free) - K(Dof_Free, Dof_Constr)*a_constr);
a(Dof_Free) = a_free;
a(Dof_Constr) = a_constr;


%% Von Mises Stresses


% Thickness Coordinate Vector
z_coord = [-h/2 0 h/2];
name = {'- h/2', '0', 'h/2'};

% Initializations
sigma_VM = zeros(n_el, length(z_coord));
figure(2)

% Cycle to loop over the different coordinates and then elements
for ii = 1:3

    % Picking z
    z = z_coord(ii);

    for el = 1:n_el

        % Displacement of the single element
        a_el = a(Edof(el, 2:end));

        % Compute Stresses
        [sigma_el, tau_el] = Stresses_Mindlin_func(Ex(el,:), Ey(el, :), a_el, z, D_Mat, G_Mat);

        % Averaging
        sigma_avg = mean(sigma_el, 2);      % Sigma 11, 22, 12
        tau_avg = mean(tau_el, 2);          % Sigma 23, 13

        % Computation of s variable to use in VM
        s = [sigma_avg(1); sigma_avg(2); 0; sigma_avg(3); tau_avg(1); tau_avg(2)] - (sigma_avg(1) + sigma_avg(2) + 0)/3*[1;1;1;0;0;0];

        % Von Mises
        sigma_VM(el, ii) = sqrt(3/2*(s'*s + sigma_avg(3)^2 + tau_avg(1)^2 + tau_avg(2)^2));
    end

    % Plotting
    subplot(1, 3, ii)
    fill(Ex', Ey', sigma_VM(:, ii)');
    title(['Von Mises Stresses at z = ', name{ii}], 'FontSize',15);
    axis equal; axis tight
    xlabel('[mm]', 'FontSize',15); ylabel('[mm]', 'FontSize',15)
    c = colorbar; c.Label.String = 'Von Mises Stress [MPa]'; c.Label.FontSize = 12;
end