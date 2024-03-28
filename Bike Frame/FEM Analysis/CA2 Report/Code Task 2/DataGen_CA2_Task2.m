function data_out = DataGen_CA2_Task2()

% Geometrical
L = 50;                          % [mm]
W = 5;                           % [mm]
h= 2.5;                          % [mm]

% Material (Linear Elastic)
E_L = 150e3;                     % [MPa]
E_T = 11e3;                      % [MPa]  
D_Mat = 1e3*[151 3.4  0;         % [MPa]
             3.4 11.1 0;
             0    0   5];


% Load
Compr_Strenght = -1.695e3;       % [MPa]
Compr_Strain = 1.13e-2;          % [1]


% Output
data_out.L = L;
data_out.W = W;
data_out.h = h;
data_out.E_T = E_T;
data_out.E_L = E_L;
data_out.D_Mat = D_Mat;
data_out.Compr_Strenght = Compr_Strenght;
data_out.Compr_Strain = Compr_Strain;

end