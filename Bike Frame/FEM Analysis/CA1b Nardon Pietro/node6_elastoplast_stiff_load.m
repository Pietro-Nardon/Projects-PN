function [K_el, fel_int] = node6_elastoplast_stiff_load(ex, ey, sigma, dsigma_deps, t)

% xi-coordinates Vector creation
xi_vect = [1/6 1/6 2/3;
           1/6 2/3 1/6];

% Integration Weights vector
H_vect = 1/6*ones(1,3);

% Creating Ke
K_el = zeros(12);
fel_int = zeros(12,1);
for ii = 1:3
    H = H_vect(ii);
    xi = xi_vect(:, ii);
    [Be, det_F_Isop] = Be_6node_func(xi, [ex(1);ey(1)], [ex(2);ey(2)], [ex(3);ey(3)], [ex(4);ey(4)], [ex(5);ey(5)], [ex(6);ey(6)]);
    K_el = K_el + Be'*dsigma_deps*Be*det_F_Isop*t*H;
    fel_int = fel_int + Be*sigma*det_F_Isop*t*H;
end





end