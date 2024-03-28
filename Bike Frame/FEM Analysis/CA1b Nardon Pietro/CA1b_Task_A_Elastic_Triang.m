%% a) - Solving plate with elasticity code
clc
clear
close all

% Datas
data_output = DataGen_CA1;

L = data_output.L;
H = data_output.H;
h = data_output.h;
E = data_output.E;
nu = data_output.poisson_ratio;
sigma_yield = data_output.sigma_yield;
H_coeff = data_output.H_coeff;

% Vector to plot to screen the results
F_React_LinearFE = zeros(5,1);

% Constitutive Matrix (Plane Strain)
D_Mat = E/((1 + nu)*(1 - 2*nu))*[1 - nu,   nu  ,      0;
                                   nu  , 1 - nu,      0;
                                   0   ,   0   , (1 - 2*nu)/2];

% Choose which mesh to load
for mesh = 1:2      % Select only the linear meshes

if mesh == 1
    load("linearCoarse2024.mat");
    name = {'Linear Coarse Mesh:   '};
elseif mesh == 2
    load("linearMedium2024.mat");
    name = {'Linear Medium Mesh:   '};
elseif mesh == 3
    load("quadraticCoarse2024.mat");
    name = {'Quadratic Coarse Mesh:'};
elseif mesh == 4
    load("quadraticMedium2024.mat");
    name = {'Quadratic Medium Mesh:'};
elseif mesh == 5
    load("quadraticFine2024.mat");
    name = {'Quadratic Fine Mesh:  '};
end


% Element nodes Coordinates
Coord_x = Coord(:,1);
Edof_x = Edof(:, 2:2:end);
Ex = Coord_x(ceil(Edof_x/2));

Coord_y = Coord(:, 2);
Edof_y = Edof(:, 3:2:end);
Ey = Coord_y(floor(Edof_y/2));


% Show Mesh
figure('WindowState','maximized')
plotpar = [1 1 0];
eldraw2_ext(Ex,Ey,plotpar); grid on
title('Structure Mesh', 'FontSize',18); subtitle(name, 'FontSize',18);


% Degrees of Freedom/Constrained (LC2)
n_elem = size(Edof, 1);
n_dof = ndofs;


% Boundary conditions
n_constr = 1 + length(bottomNodes) + length(topNodes);
Dof_Constr = [Dof(bottomNodes(1),1), reshape(Dof(bottomNodes, 2)', 1, []), reshape(Dof(topNodes, 2)', 1, [])];
Dof_Free = reshape(Dof', 1, []);
Dof_Free([Dof(bottomNodes(1),1), reshape(Dof(bottomNodes, 2)', 1, []), reshape(Dof(topNodes, 2)', 1, [])])= [];


% Initializations of Stiffness Matrix, Load Vector and Displacement Vector
K = spalloc(n_dof, n_dof, 20*n_dof);
f_ext = zeros(n_dof, 1);
a_constr = zeros(n_constr,1);    
a_constr(length(bottomNodes) + 2:end) = 1;          % Imposing a Displacement of 1mm on all the constrained nodes after the bottom ones

% Iterative cycle passing for each element
for el = 1:n_elem          
    % Computing the K and f for each element
    [K_el, f_el] = CST_stiff_load([Ex(el, 1), Ex(el,2), Ex(el,3); Ey(el,1), Ey(el,2), Ey(el,3)], D_Mat, 0, h);

    % Assemblying
    K(Edof(el, 2:end), Edof(el, 2:end)) = K(Edof(el, 2:end),Edof(el, 2:end)) + K_el;
    f_ext(Edof(el, 2:end)) = f_ext(Edof(el, 2:end)) + f_el;
end


% Solving the system equation
a_free = K(Dof_Free, Dof_Free)\(f_ext(Dof_Free) - K(Dof_Free, Dof_Constr)*a_constr);
f_ext_constr = K(Dof_Constr, Dof_Free)*a_free + K(Dof_Constr, Dof_Constr)*a_constr - f_ext(Dof_Constr);


% Total Vertical Reaction Force
F_React_LinearFE(mesh) = sum(f_ext_constr((length(bottomNodes) +2): end));
fprintf(['For the ', name{1}, ' Total Vertical Reaction Force is: ', num2str(F_React_LinearFE(mesh)), ' N\n'])


% Assemblying solutions
a(Dof_Free) = a_free;
a(Dof_Constr) = a_constr;

Q = zeros(n_dof,1);
Q(Dof_Constr) = f_ext_constr;


% Plot Deformed Shape
Ed = extract_dofs(Edof,a);         % Extract element displacements for plotting
plotpar=[2 4 0];
sfac = 1.56;                          % Magnification
eldisp2(Ex,Ey,Ed,plotpar,sfac); grid on; title('Meshed Plate'); 
xlabel('[mm]', 'FontSize', 18); ylabel('[mm]', 'FontSize', 18);axis equal

end
