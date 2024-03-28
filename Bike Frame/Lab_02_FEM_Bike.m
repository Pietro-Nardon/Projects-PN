clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bike FEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Datas

alpha = 6;                      % [1/s]
beta = 1e-5;                    % [s]
k = 4e5;                        % [N/m]
m_lumped = [1 2 2 1.5];         % [kg]
Beam_D = [20 32 50]*1e-3;       % [m]   Red - Blue - Green 
t = [1 1.5 2.5]*1e-3;           % [m]   Red - Blue - Green 
E = 7e10;                       % [N/m^2] Aluminium
rho = 2700;                     % [kg/m^3]
Omega_Max = 200*2*pi;           % [Hz] Max Freq working range

N_Springs = 2;                  % Number of Lumped Springs in the system
N_Masses = length(m_lumped);    % Number of Lumped Masses in the system


% Computations of the Sections datas
Beam_A = pi*((Beam_D).^2 - (Beam_D - 2*t).^2)/4;   % [m^2] 
Beam_J = pi*(Beam_D.^4 - (Beam_D - 2*t).^4)/64;    % [m^4]
Beam_EA = Beam_A * E;                              % [Pa*m^2]
Beam_EJ = Beam_J * E;                              % [Pa*m^4]
Beam_m = rho*(Beam_A);                             % [kg/m]


%% FEM Modelling

% Structure Assembly
[FileName, Nodes_Coord, Nodes_Num, Struct_Size, idb, N_Dof, Incid_Mat, Elem_L, Elem_Gamma, Elem_m, Elem_EA, Elem_EJ, Elem_Pos, Elem_Num, Elem_Prop] = loadstructure;
N_Dof_Total = Nodes_Num*3;

% Structure Drawing
dis_stru(Elem_Pos, Elem_L, Elem_Gamma, Nodes_Coord, Elem_Prop, idb, N_Dof);

% Matrix Assembly
[M, K] = assem(Incid_Mat, Elem_L, Elem_m, Elem_EA, Elem_EJ, Elem_Gamma, idb);


%% Lumped Parameters Implementation

% Index of the nodes that contain a Lumped Paramater
Lumped_Index_A = idb(1, 1:3);
Lumped_Index_B = idb(2, 1:3);
Lumped_Index_C = idb(3, 1:3);
Lumped_Index_F = idb(6, 1:3);

% Matrix Creation for a cycle
Lumped_Index_MAT_m = [Lumped_Index_A; Lumped_Index_B; Lumped_Index_C; Lumped_Index_F];
Lumped_Index_MAT_k = [Lumped_Index_A; Lumped_Index_B];

% Mass Matrix Computation
M_hat_MAT = zeros(3,3, N_Masses);
E_m_MAT = zeros(3, N_Dof_Total, N_Masses);
M_lumped_MAT = zeros(N_Dof_Total, N_Dof_Total, N_Masses);
for ii = 1: N_Masses
    M_hat_MAT(:, :, ii) = diag([m_lumped(ii) m_lumped(ii) 0]);
    E_m_MAT(:, Lumped_Index_MAT_m(ii,:), ii) = eye(3);                                  % Creating Expansion Matrix (3D matrix)
    M_lumped_MAT(:, :, ii) = E_m_MAT(:, :, ii)'*M_hat_MAT(:,:, ii)*E_m_MAT(:,:, ii);    % Creating Lumped Mass Matrix (3D matrix)
    M = M + M_lumped_MAT(:,:, ii);                                                      % Computing final Mass Matrix (2D matrix Dof_Tot x Dof_Tot)
end

% Stiffness Matrix Computation
K_lumped_MAT = zeros(N_Dof_Total, N_Dof_Total, N_Springs);
E_k_Vect = zeros(N_Springs, N_Dof_Total);
for ii = 1: N_Springs
    E_k_Vect(ii, Lumped_Index_MAT_k(ii, 2)) = 1;                                        % Creating Expansion Vector (2D matrix)
    K_lumped_MAT(:, :, ii) = E_k_Vect(ii, :)'*k*E_k_Vect(ii, :);                        % Creating Lumped Spring Matrix (3D matrix)
    K = K + K_lumped_MAT(:, :, ii);                                                     % Computing final Stiffness Matrix (2D matrix Dof_Tot x Dof_Tot)
end


% Damping Matrix Computation (Rayleigh Model)
C = alpha*M + beta*K;

% Free-Free Matrices extraction
M_FF = M(1:N_Dof, 1:N_Dof);
K_FF = K(1:N_Dof, 1:N_Dof);
C_FF = C(1:N_Dof, 1:N_Dof);

% Free-Contrained Matrices extraction
M_FC = M(1:N_Dof, end);
K_FC = K(1:N_Dof, end);
C_FC = C(1:N_Dof, end);

% Contrained-Free Matrices extraction
M_CF = M(end, 1:N_Dof);
K_CF = K(end, 1:N_Dof);
C_CF = C(end, 1:N_Dof);

% Constrained-Contrained Matrices extraction
M_CC = M(end, end);
K_CC = K(end, end);
C_CC = C(end, end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Requests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  1) Check for Element Lenght

% Compute L_Max for each element
L_Max_Vect = sqrt((pi^2/(2*Omega_Max))*sqrt(Elem_EJ./Elem_m));

% Print to screen L_Max check
if sum((L_Max_Vect < Elem_L)) == 0 
    fprintf('\nAll the elements work in Quasi-Static Region\n')
else 
    % Finding the position of the new nodes I need to add
    Index_Half = find(L_Max_Vect < Elem_L);
    Elem_Pos_Half = Elem_Pos + [(Elem_L.*cos(Elem_Gamma))' (Elem_L.*sin(Elem_Gamma))']/2;

    fprintf('\nThe new nodes should go in: \n')
    disp(num2str(Elem_Pos_Half(Index_Half, :)));
    
    error(['ERROR -- The elements: ' num2str(Index_Half) ' are too long'])
end



%% 2) Natural Frequencies and Mode Shapes

% Number of modes I want to analyze
N_Mode = 4;

% Find Natural Modes
[modes, omega_SQ] = eig(M_FF\K_FF);
omega = diag(sqrt(omega_SQ));                            % Bcs eig gives u the square of nat freq
[omega_sorted, omega_sorted_index] = sort((omega));      % To sort in ascending order the natural omegas
modes_sorted = modes(:, omega_sorted_index);


% Plot Mode Shapes
scale_factor = 1;
figure('WindowState','maximized')
sgtitle('            Bike Mode Shapes ', 'FontWeight', 'bold')
for ii = 1 : N_Mode
    mode = modes_sorted(:, ii);
    subplot(2,ceil(N_Mode/2), ii)
    diseg2 (mode, scale_factor, Incid_Mat, Elem_L, Elem_Gamma, Elem_Pos, idb, Nodes_Coord); 
    title(['Mode Shape: ' num2str(ii)], ['Natural Frequency: ' num2str(omega_sorted(ii)/(2*pi)) 'Hz']);
end


%% 3a) FRF Computation

% Forcing Vector
omega_F = (0: 0.1: 200)*(2*pi);     %[rad/s]  

% Find coordinates of the Input force
F0 = zeros(N_Dof, 1);
index_F0 = idb(3,2);
F0(index_F0) = 1;

% Compute FRF
X = zeros(N_Dof, length(omega_F));
for ii = 1: length(omega_F)
    A = -omega_F(ii)^2*M_FF + 1i*omega_F(ii)*C_FF + K_FF;
    X(:, ii) = A\F0;
end

% Find vertical displ and acc of F
index_F_v = idb(6,2);
FRF_F_Displ = X(index_F_v, :);
FRF_F_Acc = -(omega_F).^2.*FRF_F_Displ;

% Find horizontal displ and acc of H
index_H_h = idb(8,1); index_H_v = idb(8,2);
FRF_H_Displ = X(index_H_h, :);
FRF_H_Acc = -(omega_F).^2.*FRF_H_Displ;




%% 4) Modal Superposition approach


% Number of mode shapes considered
Mode_Shapes = 1:2;
Phi = modes_sorted(:, Mode_Shapes);


% Compute the modal matrices 
M_Modal = Phi'*M_FF*Phi;
C_Modal = Phi'*C_FF*Phi;
K_Modal = Phi'*K_FF*Phi;
F_Modal = Phi'*F0;


% FRF in modal superposition approach
XX_Modal = zeros(length(Mode_Shapes), length(omega_F));
for ii = 1 : length(omega_F)
    XX_Modal(:, ii) = (-omega_F(ii)^2 * M_Modal + 1i*omega_F(ii) * C_Modal + K_Modal) \ F_Modal;
end

% Considering just the first 2 modes
XX_MS_Modal = Phi * XX_Modal;

% Computing FRF for node F
FRF_F_Displ_Modal = XX_MS_Modal(index_F_v, :);
FRF_F_Acc_Modal = -(omega_F).^2.*FRF_F_Displ_Modal;

% Computing FRF for node H
FRF_H_Displ_Modal = XX_MS_Modal(index_H_h, :);
FRF_H_Acc_Modal = -(omega_F).^2.*FRF_H_Displ_Modal;



%% Plot of FRF numerical and Modal


FRF_Plot = [FRF_F_Displ; FRF_H_Displ; FRF_F_Acc; FRF_H_Acc];                                           % Matrix of FRF just to plot with a for cycle
FRF_Plot_Modal = [FRF_F_Displ_Modal; FRF_H_Displ_Modal; FRF_F_Acc_Modal; FRF_H_Acc_Modal];             % Matrix of FRF just to plot with a for cycle
txt_tit = {'Displacement in F', 'Displacement in H', 'Acceleration in F','Acceleration in H'};         % Titles of the plot
txt_ampl = {'[$\frac{m}{N}$]', '[$\frac{m}{N}$]','[$\frac{m}{Ns^2}$]', '[$\frac{m}{Ns^2}$]'};          % Titles of the y line of the plot
tile_plot_index = [1 2 5 6 3 4 7 8];                                                                   % Vector just to plot the graphs in a Amplitude-Phase Pattern
figure('WindowState', 'maximized');
tiledlayout(4,2);

for ii = 1 : size(FRF_Plot, 1)
    % Amplitude
    nexttile(tile_plot_index(ii))
    semilogy(omega_F/(2*pi), abs(FRF_Plot(ii, :)), 'LineWidth',2); grid on; hold on;
    semilogy(omega_F/(2*pi), abs(FRF_Plot_Modal(ii, :)), 'LineWidth',2);
    title (txt_tit(ii), 'FontSize', 15);  subtitle('Amplitude');
    xlabel('Frequency  $\Omega$ [Hz]', 'Interpreter','latex', 'FontSize', 12); ylabel (txt_ampl(ii), 'Interpreter','latex', 'FontSize',12);
    hold off;
   
    % Phase
    nexttile(tile_plot_index(ii + 4));
    plot(omega_F/(2*pi), phase(FRF_Plot(ii, :)), 'LineWidth',2); grid on; hold on
    plot(omega_F/(2*pi), phase(FRF_Plot_Modal(ii, :)), 'LineWidth',2); 
    title ('', 'Phase'); 
    xlabel('Frequency  $\Omega$ [Hz]', 'Interpreter','latex', 'FontSize', 12); ylabel('$\phi$ [rad]', 'Interpreter','latex', 'FontSize',12);
    yline(0, 'k--');
    hold off;
end

% Legend of the Plots
Lgnd = legend('FEM Approach', 'Modal Approach');
Lgnd.Position(1) = 0.465;
Lgnd.Position(2) = 0.94;



%% 3b) Internal Forces of midpoint of GE

% Numbers of the beams that touch the midpoint
n_el_dx = 4;
n_el_sx = 5;
n_el = [n_el_dx n_el_sx];
L_el = Elem_L(n_el);

% Number of the nodes
index_dof_i = [idb(9,:); idb(7,:)];
index_dof_j = [idb(5,:); idb(9,:)];

% Axial Position
ax_pos = [0 L_el(2)];

N_IF = zeros(length(n_el), size(X, 2));
T_IF = zeros(length(n_el), size(X, 2));
M_IF = zeros(length(n_el), size(X, 2));
for ii = 1:2
    % Rotational Matrix
    Rot_Matrix = [cos(Elem_Gamma(n_el(ii))) sin(Elem_Gamma(n_el(ii))) 0;
                 -sin(Elem_Gamma(n_el(ii))) cos(Elem_Gamma(n_el(ii))) 0; 
                         0                    0             1];

    % Rot_Matrix = eye(3);
    
    % Computing FRF local node
    Xi = Rot_Matrix*X(index_dof_i(ii, :), :);
    Xj = Rot_Matrix*X(index_dof_j(ii, :), :);
    
    % Computing Coefficients
    b = (Xj(1,:) - Xi(1,:))/L_el(ii);
    c = -3/L_el(ii)^2*Xi(2,:) + 3/L_el(ii)^2*Xj(2, :) - 2/L_el(ii)*Xi(3, :) - 1/L_el(ii)*Xj(3, :);
    d =  2/L_el(ii)^3*Xi(2,:) - 2/L_el(ii)^3*Xj(2, :) + 1/L_el(ii)^2*Xi(3, :) + 1/L_el(ii)^2*Xj(3, :);

    % Computing Internal Forces
    N_IF(ii, :) = Elem_EA(n_el(ii))*b;
    T_IF(ii, :) = Elem_EJ(n_el(ii))*(6*d);
    M_IF(ii, :) = Elem_EJ(n_el(ii))*(2*c + 6*d*ax_pos(ii));
end

Int_Forces = [N_IF;T_IF;M_IF];

% Plots
txt_tit = {'Axial Force N','','Shear Force T','', 'Bending Moment M'};
txt_ampl = {'$|\frac{N}{F}|$ [$\frac{N}{N}$]', '', '$|\frac{T}{F}|$ [$\frac{N}{N}$]', '', '$|\frac{M}{F}|$ [$\frac{N}{N}$]', ''};

tile_plot_index = [1 4 2 5 3 6];
figure('WindowState', 'maximized'); tiledlayout(2,3);
for ii = [1 3 5]
    nexttile(tile_plot_index(ii))
    semilogy(omega_F/(2*pi), abs(Int_Forces(ii, :)), omega_F/(2*pi), abs(Int_Forces(ii + 1, :)), 'LineWidth',2); grid on;
    title(txt_tit(ii)); subtitle('Amplitude'); legend('Right Element', 'Left Element')
    xlabel('Frequency  $\Omega$ [Hz]', 'Interpreter','latex', 'FontSize', 12); ylabel(txt_ampl(ii), 'Interpreter','latex', 'FontSize',12);

    nexttile(tile_plot_index(ii + 1))
    plot(omega_F/(2*pi), phase(Int_Forces(ii, :)), omega_F/(2*pi), phase(Int_Forces(ii + 1, :)), 'LineWidth',2); grid on;
    subtitle('Phase')
    xlabel('Frequency  $\Omega$ [Hz]', 'Interpreter','latex', 'FontSize', 12); ylabel('$\phi$ [rad]', 'Interpreter','latex', 'FontSize',12); 
    yline(0, 'k--'); ylim([-5 3])
end
sgtitle('Internal Forces in midpoint of GE', 'FontSize', 15, 'FontWeight', 'bold')

%% 3c) Contrain Forces
% all the xc will be zero while the Ff = F0 and Fc = 0
R = zeros(1, length(omega_F)); X_Constr = zeros(N_Dof, length(omega_F));
for ii = 1:length(omega_F)
    % X_Constr(:, ii) = -(-M_FF*omega_F(ii)^2 + 1i*C_FF*omega_F(ii) + K_FF);   No need bcs it's the same as computing X 
    R(:, ii) = (-M_CF*omega_F(ii)^2 + 1i*C_CF*omega_F(ii) + K_CF)*X(:, ii);
end

% Plot
figure('WindowState', 'maximized'); 
subplot(211)
semilogy(omega_F/(2*pi), abs(R), 'LineWidth',2); grid on;
xlabel('Frequency  $\Omega$ [Hz]', 'Interpreter','latex', 'FontSize', 12); ylabel('$|R|$ $[N]$', 'Interpreter','latex', 'FontSize',12); 
title('Constraint Force in C', 'Amplitude');

subplot(212)
plot(omega_F/(2*pi), phase(R), 'LineWidth',2); grid on;
xlabel('Frequency  $\Omega$ [Hz]', 'Interpreter','latex', 'FontSize', 12); ylabel('$\phi$  [rad]', 'Interpreter','latex', 'FontSize', 12);
subtitle('Phase'); 
yline(0, 'k--')



%% 5) FRF of F with Input Vertical Displ in A


F_Spring = zeros(N_Dof, 1);
index_A_v = idb(1,2);
F_Spring(index_A_v) = k;            % I compute the moving of A' by just finding the contribution of the movement of the bottom end
                                    % it should be k*y but being a FRF I just put y = 1
                               

% Compute FRF
X_Spring = zeros(N_Dof, length(omega_F));
for ii = 1: length(omega_F)
    X_Spring(:, ii) = (-omega_F(ii)^2*M_FF + 1i*omega_F(ii)*C_FF + K_FF)\F_Spring;
end

% Find the FRF of (F,2)
FRF_Spring = -(omega_F.^2).*X_Spring(index_F_v, :);

% Plot
figure('WindowState', 'maximized'); 
subplot(211)
semilogy(omega_F/(2*pi), abs(FRF_Spring), 'LineWidth',2); grid on;
xlabel('Frequency  $\Omega$ [Hz]', 'Interpreter','latex', 'FontSize', 12); ylabel('Acceleration [$\frac{m}{Ns^2}$]', 'Interpreter','latex', 'FontSize', 12); %xlim([1 200]) % to make the graph prettier
title('Vertical Acceleration point F', 'Amplitude');

subplot(212)
plot(omega_F/(2*pi), phase(FRF_Spring), 'LineWidth',2); grid on;
xlabel('Frequency  $\Omega$ [Hz]', 'Interpreter','latex', 'FontSize', 12); ylabel('$\phi$  [rad]', 'Interpreter','latex', 'FontSize', 12);
subtitle('Phase'); 
yline(0, 'k--');




%% 6) Time history irregular surface

% Datas
bike_speed = 12;    %[m/s]
lambda_1 = 1;       %[m]
lambda_2 = 0.6;     %[m]
A_1 = 1e-3;         %[m]
A_2 = 5e-4;         %[m]
phi_1 = 0;          %[rad]
phi_2 = 0;          %[rad]

% Forcing frequency
f_1 = bike_speed/lambda_1;
f_2 = bike_speed/lambda_2;

% Angular Speeds
omega_1 = 2*pi*f_1;
omega_2 = 2*pi*f_2;

% Periods 
T_1 = 1/f_1;
T_2 = 1/f_2;

% Bike time
bike_pitch = Nodes_Coord(1,1) - Nodes_Coord(2,1);       %[m]
bike_time = bike_pitch/bike_speed;                      %[s]

% Phase delay
delay_1 = omega_1*(bike_time - floor(bike_time/T_1)*T_1);      % Floor to have a single phase delay, not an entire circle
delay_2 = omega_2*(bike_time - floor(bike_time/T_2)*T_2);      % Floor to have a single phase delay, not an entire circle


% Forces Applied
F_Road_1 = zeros(N_Dof, 1);
F_Road_2 = zeros(N_Dof, 1);
index_B_v = idb(2,2);

F_Road_1 ([index_A_v, index_B_v]) = [A_1*k, A_1*k*exp(1i*delay_1)];
F_Road_2 ([index_A_v, index_B_v]) = [A_2*k, A_2*k*exp(1i*delay_2)];


% Compute FRF
FRF_Road_1 = zeros(N_Dof, length(omega_F));
FRF_Road_2 = zeros(N_Dof, length(omega_F));
for ii = 1: length(omega_F)
    FRF_Road_1(:, ii) = (-omega_F(ii)^2*M_FF + 1i*omega_F(ii)*C_FF + K_FF)\F_Road_1;
    FRF_Road_2(:, ii) = (-omega_F(ii)^2*M_FF + 1i*omega_F(ii)*C_FF + K_FF)\F_Road_2;
end

% Compute Acceleration
FRF_Road_Acc_1_H = abs(-(omega_F.^2).*FRF_Road_1(index_H_v, :));
FRF_Road_Acc_2_H = abs(-(omega_F.^2).*FRF_Road_2(index_H_v, :));

% Amplitude of the FRF at their omega
Amplitude_1 = FRF_Road_Acc_1_H(omega_1/(0.1*2*pi));
Amplitude_2 = FRF_Road_Acc_2_H(omega_2/(0.1*2*pi));

% Time Response
t_vector = 0:0.001:1;   %[s]
time_response = Amplitude_1*cos(omega_1*t_vector) + Amplitude_2*cos(omega_2*t_vector);

% Plot
figure('WindowState','maximized')
plot(t_vector, time_response, 'LineWidth',2); grid on;
title('Vertical Acceleration of Point H due to road irregularity')
xlabel('Time $[s]$', 'interpreter', 'latex', 'FontSize',12); ylabel('Acceleration [$\frac{m}{s^2}$]', 'Interpreter','latex', 'FontSize',12);
yline(0, 'k--')




%% 7) Static Displacement

% Creation of the Static Force
F_Static = zeros(N_Dof, 1);
F_Static(index_H_v) = -600;      % Saddle Force [N]
F_Static(index_F_v) = -100;      % Handle Force [N]

% Computic static displacement
Static_Displ = K_FF \ F_Static(1:N_Dof);

% Plotting Static Deformation
figure('WindowState','maximized')
diseg2(Static_Displ, 100, Incid_Mat, Elem_L, Elem_Gamma, Elem_Pos, idb, Nodes_Coord);
title('Bike Static Deformation')

% Compute vertical displacement
index_9_v = idb(9,2);
Vert_Displ_9 = Static_Displ(index_9_v);
fprintf(['\nThe vertical displacement of the midpoint of GE is: ' num2str(Vert_Displ_9) ' mm\n'])






% close all





