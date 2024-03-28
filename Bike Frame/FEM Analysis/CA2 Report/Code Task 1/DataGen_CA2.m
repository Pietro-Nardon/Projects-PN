function data_out = DataGen_CA2()

% Geometrical
L = 30e1;             % [mm]
W = 20e1;             % [mm]
h_min = 0.1e1;        % [mm]
h_max = 10e1;         % [mm]

% Material (Linear Elastic)
E_Mat = 3e3;          % [MPa]
nu = 0.3;             % [1]  

% Load
p = -10e-1;           % [MPa]



% Output
data_out.L = L;
data_out.W = W;
data_out.h_min = h_min;
data_out.h_max = h_max;
data_out.E_Mat = E_Mat;
data_out.nu = nu;
data_out.p = p;

end