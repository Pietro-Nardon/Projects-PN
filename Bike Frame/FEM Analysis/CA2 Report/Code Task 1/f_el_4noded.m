function [f_el_ext] = f_el_4noded(Ex, Ey, p)
% This function computes the element vertical load contribution due to
% pressure

% xi-coordinates Vector and Weights
H_vect = ones(1,4);
xi_vect = [-1/sqrt(3) -1/sqrt(3)  1/sqrt(3)  1/sqrt(3);
           -1/sqrt(3)  1/sqrt(3) -1/sqrt(3)  1/sqrt(3)];

% Initializing
f_el_ext = zeros(4,1);


% Loop over integration points
for ii = 1:4 
    H = H_vect(ii);
    xi = xi_vect(:, ii); 
    [~, det_Fisop, N_w] = Vert_Be_quad_func(xi, [Ex(1);Ey(1)], [Ex(2);Ey(2)], [Ex(3);Ey(3)], [Ex(4);Ey(4)]);
    f_el_ext = f_el_ext + N_w'*p*det_Fisop*H;
end
end