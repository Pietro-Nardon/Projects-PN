%% d) - ElastoPlastic Behaviour
clc
clear
close all

tic
% Datas
data_output = DataGen_CA1;

L = data_output.L;
H = data_output.H;
h = data_output.h;
E = data_output.E;
nu = data_output.poisson_ratio;
sigma_yield = data_output.sigma_yield;
H_coeff = data_output.H_coeff;

F_React_NonLinearFE = zeros(5, 1);
num_iter = zeros(100,3);


% Input vectors for CALFEM (They won't change in the loop)
ptype = 2;
Mat_Prop = [E, nu, H_coeff];
D_Mat = hooke(ptype, E, nu);

% Choose which mesh to load
for mesh = 3:5


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



    % Degrees of Freedom/Constrained (LC2)
    n_elem = size(Edof, 1);
    n_dof = ndofs;
    n_constr = 1 + length(bottomNodes) + length(topNodes);
    Dof_Prescribed = reshape(Dof(topNodes, 2)', 1, []);
    Dof_Constr = [Dof(bottomNodes(1),1), reshape(Dof(bottomNodes, 2)', 1, []), reshape(Dof(topNodes, 2)', 1, [])];
    Dof_Free = reshape(Dof', 1, []);
    Dof_Free([Dof(bottomNodes(1),1), reshape(Dof(bottomNodes, 2)', 1, []), reshape(Dof(topNodes, 2)', 1, [])])= [];


    % Initialize Displacements
    a = zeros(n_dof, 1);
    a_old = zeros(n_dof, 1);
    d_a = a - a_old;
    a_constr = zeros(n_constr,1);
    
    % Initialize State Variables
    n_state = 3;
    ngp = 3;                            % Number of Gauss integration points    
    stress_old = zeros(4, ngp, n_elem);
    stress = stress_old;
    state = zeros(n_state, ngp, n_elem);
    state_old = zeros(size(state));

    % Imposing yield stress
    state_old(2,:,:) = sigma_yield;

    % Time Stepping
    n_time = 250;               % Number of timesteps
    t_end = 100;                
    t = linspace(0, t_end, n_time);

    % Displacement Control
    U_max = 1;                  %[mm]
    Uu = linspace(0, U_max, n_time);

    % Initialize Variables for post processing
    K = spalloc(n_dof, n_dof, 20*n_dof);
    F_react = zeros(size(Uu));

    % Inizialize force vectors
    f_int = zeros(n_dof, 1);
    f_ext = f_int;

    % Gravity Load
    F_Body = zeros(n_elem, 2);    % In this case is absent

    % Tolerance for Newton iterations
    toll = 1e-6;

    figure('WindowState','maximized')
    for ii = 1: n_time

        % Initial Guess of unknown Displecement field
        a(Dof_Free) = a_old(Dof_Free) + d_a(Dof_Free);

        % Update Prescribed dof
        a(Dof_Prescribed) = Uu(ii);

        % Iterations with Newton
        unbal = 1e10;         % It's the error 
        n_iter = 0;
        stress_trial = zeros(4,1);

        while unbal > toll
            % Nullifying
            K = K.*0;
            f_int = f_int.*0;
        
            % Loop over the elements
            for el = 1: n_elem
                % Find dof of the element
                dof = Edof(el, 2:end);

                % Compute incrementral element displacement
                d_ae = a(dof) - a_old(dof);

                % Compute Stiffness Mat, Force vect, stress and state vect
                [K_el, f_el, stress(:,:, el), state(:,:, el)] = node6_isop_non_lin(Ex(el,:), Ey(el,:), d_ae, h, D_Mat, ptype,Mat_Prop,stress_old(:,:,el),state_old(:,:,el));
               
                % Assemblying 
                K(Edof(el, 2:end), Edof(el, 2:end)) = K(Edof(el, 2:end), Edof(el, 2:end)) + K_el;
                f_int(Edof(el, 2:end)) = f_int(Edof(el, 2:end)) + f_el;
              
            end
            % Compute Unbalance
            g_F = f_int(Dof_Free) - f_ext(Dof_Free);
            unbal = norm(g_F);
    
            % Newton Update
            if unbal > toll
                a(Dof_Free) = a(Dof_Free) - K(Dof_Free, Dof_Free)\g_F;
            end

            % Update Iterations
            n_iter = n_iter + 1;
    
            % Divergence Check
            if n_iter > 20
                error('Newton cannot converge')
            end

        end

        % Computing external Forces
        f_ext_constr = f_int(Dof_Constr);


        % Save Data on Reaction Force
        F_react(ii) = sum(f_ext_constr((length(bottomNodes) + 2): end));

        % Save New displacement
        d_a = a - a_old;
        a_old = a;

        % Storing
        num_iter(ii, mesh - 2) = n_iter;
        state_old = state;
        stress_old = stress;

       
        % % Plots
        % plot(Uu, F_react, '-', 'LineWidth',2); grid on; 
        % xlabel('Load Time Stepping [s]'); ylabel('Reaction Force [N]');
        % title('Vertical Reaction developement over Prescribed Displacement')
        % subtitle(name); ylim([0 4e5]);
        % drawnow
    end
    % Screen Display
    fprintf(['For the ', name{1}, ' Total Vertical Reaction Force is: ', num2str(F_react(end)), ' N\n'])
    F_React_NonLinearFE(mesh) = F_react(end);

    % Plots
    plot(Uu, F_react, '-', 'LineWidth',2); grid on; 
    xlabel('Load Time Stepping [s]'); ylabel('Reaction Force [N]');
    title('Vertical Reaction developement over Prescribed Displacement')
    subtitle(name); ylim([0 4e5]);

    % Plots
    figure('WindowState','maximized')
    plotpar = [1 1 0];
    eldraw2_ext(Ex,Ey,plotpar); grid on
    xlabel('[mm]'); ylabel('[mm]'); 
    title('Structure Mesh'); subtitle(name);
    Ed = extract_dofs(Edof,a);
    plotpar=[2 4 0];
    sfac = 1;                         
    eldisp2_ext(Ex,Ey,Ed,plotpar,sfac); grid on;


% End of the for loop to do more mesh at same time 
end
toc
n_time