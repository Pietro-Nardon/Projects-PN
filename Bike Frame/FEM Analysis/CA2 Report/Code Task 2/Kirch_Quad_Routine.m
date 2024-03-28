function [K_uu_el, K_ww_K_el, G_el, fe_w] = Kirch_Quad_Routine(Ex, Ey, D_Mat, N_sec, h, p)
% Output:    K_uu_el =  [8x8]
%            K_ww_el = [12x12]               
%               G_el = [12x12]  
%               fe_u = [12x1]


% Material Matrices
D_bar = h^3/12*D_Mat;
D_0 = h*D_Mat;

%% Cycle over 4 Integration Points

% Gauss Integration points and weihts
H_vect_4 = ones(1,4);
xi_vect_4 = [-1/sqrt(3) -1/sqrt(3)  1/sqrt(3)  1/sqrt(3);
             -1/sqrt(3)  1/sqrt(3) -1/sqrt(3)  1/sqrt(3)];

% Initializations
K_uu_el = zeros(8);
K_ww_K_el = zeros(12);
fe_w=zeros(8,1);
P = [0; p];

for ii = 1:4
    H = H_vect_4(ii);
    xi = xi_vect_4(:,ii);

    % Computation of determinant and B matrices
    det_Fisop = detFisop_4node_func(xi,[ Ex(1) Ey(1) ]',[ Ex(2) Ey(2) ]',[ Ex(3) Ey(3) ]',[ Ex(4) Ey(4) ]');
    Bastn = Bast_kirchoff_func(xi,[ Ex(1) Ey(1) ]',[ Ex(2) Ey(2) ]',[ Ex(3) Ey(3) ]',[ Ex(4) Ey(4) ]');
    Be_u = Be_Mindlin_func(xi,[ Ex(1) Ey(1) ]',[ Ex(2) Ey(2) ]',[ Ex(3) Ey(3) ]',[ Ex(4) Ey(4) ]');
    Ne_u = N_Mindlin_func(xi,[ Ex(1) Ey(1) ]',[ Ex(2) Ey(2) ]',[ Ex(3) Ey(3) ]',[ Ex(4) Ey(4) ]');

    % Computation of Matrices
    K_uu_el = K_uu_el + Be_u'*D_0*Be_u*det_Fisop*H;
    K_ww_K_el = K_ww_K_el + Bastn'*D_bar*Bastn*det_Fisop*H;

    % Computation of the force
    fe_w = fe_w + Ne_u'*P*det_Fisop*H;

end


%% Cycle over 9 integration points

% Gauss Integration points and weights
H_vect_9 = [0.309 0.494 0.309 0.494 0.790 0.494 0.309 0.494 0.309];
xi_vect_9 = [-sqrt(3/5)      0      sqrt(3/5) -sqrt(3/5)    0    sqrt(3/5) -sqrt(3/5)     0     sqrt(3/5);
             -sqrt(3/5) -sqrt(3/5) -sqrt(3/5)      0        0        0      sqrt(3/5) sqrt(3/5) sqrt(3/5)];

% Initialization
G_el = zeros(12);

for ii = 1:9
    H = H_vect_9(ii);
    xi = xi_vect_9(:,ii);

    % Compute determinant B_el Matrix
    det_Fisop = detFisop_4node_func(xi,[ Ex(1) Ey(1) ]',[ Ex(2) Ey(2) ]',[ Ex(3) Ey(3) ]',[ Ex(4) Ey(4) ]');
    [~, B_el_T] = Bast_kirchoff_func(xi,[ Ex(1) Ey(1) ]',[ Ex(2) Ey(2) ]',[ Ex(3) Ey(3) ]',[ Ex(4) Ey(4) ]');
    
    % Compute G matrix
    G_el = G_el + B_el_T*N_sec*B_el_T'*det_Fisop*H;
end















end
