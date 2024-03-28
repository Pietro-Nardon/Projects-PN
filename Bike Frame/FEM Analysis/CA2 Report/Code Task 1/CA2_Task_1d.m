clc
clear
close all

%% %%%%%%%%%%%%%%% Computer Assignment 2 - Task_1 d) %%%%%%%%%%%%%%%%%%%%%%


%% Loading Datas
data = DataGen_CA2;

L = data.L;
W = data.W;
h_min = data.h_min;
h_max = data.h_max;
E_Mat = data.E_Mat;
nu = data.nu;
p = data.p;

% Costitutive Matrix
ptype = 1;              % Plane Stress
D_Mat = hooke(ptype, E_Mat, nu);

%% -------- Kirchhoff Plate --------------------------------------

% Creating the mesh
n_el_x = 20;
n_el_y = 30;
n_type = 3;             % 3 Kirch       % 5 Mindlin
[mesh, Coord, Edof]=rectMesh(0, W, 0, L, n_el_x, n_el_y, n_type);

% Element nodes Coordinates
Coord_x = Coord(:,1);
Edof_x = Edof(:, 2:n_type:end);
Ex = Coord_x(ceil(Edof_x/n_type));

Coord_y = Coord(:, 2);
Edof_y = Edof(:, 3:n_type:end);
Ey = Coord_y(ceil(Edof_y/n_type));

% Plot mesh
figure(1)
subplot(121)
plotpar = [1 1 0];
eldraw2(Ex,Ey,plotpar); grid on
title('Undeformed Plate with Kirchhoff Elements', 'FontSize',15); axis tight
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

Dof = reshape(1:n_dof, n_type, [])';

% Applying Contrains 
Dof_Free = 1:n_dof;
Dof_Constr_w = [Dof(Nodes_Bottom, 1); Dof(Nodes_Right,1); Dof(Nodes_Left,1); Dof(Nodes_Top,1)];
Dof_Constr_Rot = [Dof(Nodes_Bottom, 2); Dof(Nodes_Right,3); Dof(Nodes_Left,3); Dof(Nodes_Top,2)];
Dof_Constr = unique([Dof_Constr_w; Dof_Constr_Rot]);
n_constr = length(Dof_Constr);
Dof_Free(Dof_Constr) = [];




%% Displacements

% Thickness Vector
h_vect = h_min:1:h_max;

% Initializations
Max_Vert_Disp_kirch_quantity = zeros(1, length(h_vect));

% Cycle to pass by every thickness
for h = h_vect

    % D_bar computation
    D_bar = h^3/12*D_Mat;

    % Resetting the variables
    clear K_el f_el_ext
    K = spalloc(n_dof,n_dof,20*n_dof);
    f_ext = zeros(n_dof,1);
    a = zeros(n_dof, 1);
    a_constr = zeros(n_constr, 1);

    % Cycle over the elements
    for el = 1:n_el
        [K_el, f_el_ext] = Kirch_K_f_4noded(Ex(el,:), Ey(el, :), p, D_bar);
        [K, f_ext] = assem(Edof(el, :), K, K_el, f_ext, f_el_ext);
    end

    % Solving system of equations
    a_free = K(Dof_Free, Dof_Free)\(f_ext(Dof_Free) - K(Dof_Free, Dof_Constr)*a_constr);
    a(Dof_Free) = a_free;
    a(Dof_Constr) = a_constr;

    % Compute and store maximum displacement
    Max_Vert_Disp_kirch_quantity(h) = min(a(1:3:end))*h^3;


end

% Plotting Kirchhoff displacements
figure(2)
subplot(223)
plot(h_vect, Max_Vert_Disp_kirch_quantity); grid on
ylim([-9 -4]*1e4); xlim tight
title('Plot of the quantity wh^3 for a Kirchhoff element plate', 'FontSize', 12);
ylabel('wh^3 [mm^4]', 'FontSize', 12); xlabel('h [mm]' , 'FontSize', 12);
subplot(221)
Ed = extract_dofs(Edof,a); % extract element displacements for plotting
fill(Ex',Ey',Ed(:,1:3:end)')
grid on; colorbar; axis equal; axis tight
title('Kirchhoff Plate Vertical Displacement [m]', 'FontSize', 15)
xlabel('[mm]', 'FontSize', 15); ylabel('[mm]', 'FontSize', 15)

figure(3)
plot(h_vect, Max_Vert_Disp_kirch_quantity, 'LineWidth',2); grid on; hold on
ylim([-9 -4]*1e4); xlim tight


%% ---------------- Mindlin ----------------------------------------

% Creating the mesh
n_el_x = 20;
n_el_y = 30;
n_type = 5;             % 3 Kirch       % 5 Mindlin
[mesh, Coord, Edof]=rectMesh(0, W, 0, L, n_el_x, n_el_y, n_type);

% Element nodes Coordinates
Coord_x = Coord(:,1);
Edof_x = Edof(:, 2:n_type:end);
Ex = Coord_x(ceil(Edof_x/n_type));

Coord_y = Coord(:, 2);
Edof_y = Edof(:, 3:n_type:end);
Ey = Coord_y(ceil(Edof_y/n_type));

% Plot mesh
figure(1)
subplot(122)
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

Dof = reshape(1:n_dof, n_type, [])';

% Applying Contrains 
Dof_Free = 1:n_dof;
Dof_Constr_x = [Dof(Nodes_Bottom, 1); Dof(Nodes_Right,1); Dof(Nodes_Left,1); Dof(Nodes_Top,1)];
Dof_Constr_y = [Dof(Nodes_Bottom, 2); Dof(Nodes_Right,2); Dof(Nodes_Left,2); Dof(Nodes_Top,2)];
Dof_Constr_w = [Dof(Nodes_Bottom, 3); Dof(Nodes_Right,3); Dof(Nodes_Left,3); Dof(Nodes_Top,3)];
Dof_Constr_Rot = [Dof(Nodes_Bottom, 4); Dof(Nodes_Right,5); Dof(Nodes_Left,5); Dof(Nodes_Top,4)];
Dof_Constr = unique([Dof_Constr_x; Dof_Constr_y; Dof_Constr_w; Dof_Constr_Rot]);
n_constr = length(Dof_Constr);
Dof_Free(Dof_Constr) = [];




%% Displacements

% Thickness Vector
h_vect = h_min:1:h_max;

% Initializations
Max_Vert_Disp_mindlin_quantity = zeros(1, length(h_vect));

% Cycle to pass by every thickness
for h = h_vect

    % Material Matrices computations
    D_bar = h^3/12*D_Mat;
    D_0 = h*D_Mat;
    G_Mat = E_Mat/(2*(1+nu));
    G_0 = G_Mat*5/6*h;

    % Resetting the variables
    clear K_el f_el_ext
    K = spalloc(n_dof,n_dof,20*n_dof);
    f_ext = zeros(n_dof,1);
    a = zeros(n_dof, 1);
    a_constr = zeros(n_constr, 1);

    % Cycle over the elements
    for el = 1:n_el
        [K_el, ~, ~, ~, ~, ~, f_el_ext] = Kel_fel_Mindlin_func(Ex(el,:), Ey(el, :), D_bar, D_0, G_0, p);
        [K, f_ext] = assem(Edof(el, :), K, K_el, f_ext, f_el_ext);
    end

    % Solving system of equations
    a_free = K(Dof_Free, Dof_Free)\(f_ext(Dof_Free) - K(Dof_Free, Dof_Constr)*a_constr);
    a(Dof_Free) = a_free;
    a(Dof_Constr) = a_constr;

    % Compute and store maximum displacement
    Max_Vert_Disp_mindlin_quantity(h) = min(a(3:5:end))*h^3;
end

% Plotting Mindlin displacements
figure(2)
subplot(224)
plot(h_vect, Max_Vert_Disp_mindlin_quantity); grid on; xlim tight
title('Plot of the quantity wh^3 for a Mindlin element plate', 'FontSize', 12);
ylabel('wh^3 [mm^4]', 'FontSize', 12); xlabel('h [mm]' , 'FontSize', 12);
subplot(222)
Ed = extract_dofs(Edof,a); 
fill(Ex',Ey',Ed(:,3:5:end)')
grid on; colorbar; axis equal; axis tight
title('Mindlin Plate Vertical Displacement [m]', 'FontSize', 15)
xlabel('[mm]', 'FontSize', 15); ylabel('[mm]', 'FontSize', 15)

figure(3)
plot(h_vect, Max_Vert_Disp_mindlin_quantity, 'LineWidth',2); grid on; xlim tight
title('Plot of the quantity w_{max}h^3', 'FontSize', 12);
ylabel('w_{max}h^3 [mm^4]', 'FontSize', 12); xlabel('h [mm]' , 'FontSize', 12);
legend('Kirchhoff', 'Mindlin', 'Location','southwest');


