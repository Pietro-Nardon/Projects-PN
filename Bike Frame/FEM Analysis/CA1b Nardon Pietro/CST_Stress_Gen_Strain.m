function Sigma_Princ = CST_Stress_Gen_Strain(a_vect, node_coord, E, nu)
%
% Calculates the Principal Stresses in the case of a CST with plane strain
% -----------------------------------------------------------------------
% Input:    
%                 a_vect = Vector of the displacements [6x1]
%             node_coord = 1 row = x; 2 row = y        [2x3]
%                      E = Elastic Module of the material
%                     nu = Poisson's Ratio of the material
% -----------------------------------------------------------------------
% Output: 
%            Sigma_Princ = [3x1]

% Costitutive Matrix (Plane Strain)
D_Mat = E/((1 + nu)*(1 - 2*nu))*[1 - nu,   nu  ,      0;
                                   nu  , 1 - nu,      0;
                                   0   ,   0   , (1 - 2*nu)/2];

% Find Be
Be = Be_cst_func(node_coord(:, 1), node_coord(:, 2),node_coord(:, 3));

% Compute Strains
Et = Be*a_vect;

% Compute Stresses
Es = D_Mat*Be*a_vect;

% Compute Sigma_zz
Sigma_zz = (nu*E)*(Et(1) + Et(2))/((1 + nu)*(1 - 2*nu));

% Create Sigma Matrix
Sigma_Mat = [Es(1) Es(3) 0;
             Es(3) Es(2) 0;
              0     0  Sigma_zz];

% Find Principal Stresses
Sigma_Princ = sort(eig(Sigma_Mat), 'descend')';



end