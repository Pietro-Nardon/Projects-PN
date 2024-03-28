function data_output = DataGen_CA1()

% Geometrical
L = 100;               % [mm]
H = 100;               % [mm]
h = 10;                % [mm]

% Material
E = 210e3;              % [MPa]
poisson_ratio = 0.3;    % [1]
sigma_yield = 500;      % [MPa]
H_coeff = 20000;        % [MPa]


% Output
data_output.L = L;
data_output.H = H;
data_output.h = h;
data_output.E = E;
data_output.poisson_ratio = poisson_ratio;
data_output.sigma_yield = sigma_yield;
data_output.H_coeff = H_coeff;

end