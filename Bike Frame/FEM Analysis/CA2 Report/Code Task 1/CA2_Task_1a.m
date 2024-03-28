clc
clear 
close all

%% %%%%%%%%%%%%%%% Computer Assignment 2 - Task_1 a) %%%%%%%%%%%%%%%%%%%%%%


%% Loading Datas
data = DataGen_CA2;

L = data.L;
W = data.W;
p = data.p;


%% Task_1 - a) Element Load Contribution

% Creating the mesh
n_el_x = 20;
n_el_y = 30;
n_type = 5;             % 5 Mindlin
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
title('Undeformed Plate', 'FontSize',15); axis tight
xlabel('[mm]', 'FontSize',15); ylabel('[mm]', 'FontSize',15)

% Numbers
n_el = size(mesh, 2);
n_nodes = size(Coord, 1);

% Loop over the elements to compute the element force due to the pressure
f_el_ext_w = zeros(4, n_el);
for el = 1: n_el
    f_el_ext_w(:, el) = f_el_4noded(Ex(el,:), Ey(el, :), p);
end

% Checking the total pressure
F_Load_sum = sum(sum(f_el_ext_w));
F_Load = p*L*W;
if (abs(F_Load - F_Load_sum) < 1e-6)
    fprintf(['\nThe sum of all load contributions is equal to the total force and it is: ', num2str(F_Load_sum), ' N\n'])

end
