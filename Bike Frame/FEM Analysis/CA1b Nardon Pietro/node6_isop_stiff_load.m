function [K_el, fel_ext] = node6_isop_stiff_load(ex, ey, D_Mat, F_Body, t)

% xi-coordinates Vector creation
xi_vect = [1/6 1/6 2/3;
           1/6 2/3 1/6];

% Integration Weights vector
H_vect = 1/6*ones(1,3);

% Initializing 
K_el = zeros(12);
fel_ext = zeros(12,1);

% Loop over integration points
for ii = 1:3
    H = H_vect(ii);
    xi = xi_vect(:, ii);
    [Be, det_F_Isop, Ne] = Be_6node_func(xi, [ex(1);ey(1)], [ex(2);ey(2)], [ex(3);ey(3)], [ex(4);ey(4)], [ex(5);ey(5)], [ex(6);ey(6)]);
    K_el = K_el + Be'*D_Mat*Be*det_F_Isop*t*H;
    fel_ext = fel_ext + Ne'*F_Body*det_F_Isop*t*H;
    
end





end