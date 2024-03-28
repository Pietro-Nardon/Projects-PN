function [K_el,f_el_int,stress,state] = node6_isop_non_lin(ex, ey, d_ae, t, D_Mat, ptype,Mat_Prop,stress_old,state_old)
% =========================================================================
% CALL: [Ke,fe_int,stress,state] = own_tri_elfun(ex, ey, d_ae, t,...
%                                       D, ptype,mp,stress_old,state_old)
%
% PURPOSE: Computes element stiffness matrix, external force vector as well
% as stress and state variables for a 6-noded isoparametric
% triangular element. von Mises plasticity with linear isotropic
% hardening is assumed.
%
% INPUT:    ex,ey       - Element coordinates.
%           d_ae        - Incremental element displacement
%           D_Mat           - Constitutive elastic matrix.
%           t           - Thickness.
%           ptype       - 1: plane stress, 2: plane strain.
%           Mat_Prop    - Material properties (E,v,H).
%                         E: modulus of elasticity
%                         v: Poisson's ratio
%                         H: modulus of plasticity
%           stress_old  - Stress from previous time step.
%           state_old   - state variables from previous time step.
%
% OUTPUT:   Ke          - Element stiffness matrix.
%           fe_int      - Element internal force vector.
%           stress      - Updated stress state.
%           state       - Updated state. 
% =========================================================================


% Define integration point and weight
ngp = 3;

% xi-coordinates Vector creation
xi_vect = [1/6 1/6 2/3;
           1/6 2/3 1/6];

% Integration Weights vector
H_vect = 1/6*ones(1,3);

% Initializing 
K_el = zeros(12,12);
f_el_int = zeros(12, 1);
state = zeros(size(state_old));
stress = zeros(size(stress_old));

% Loop over Gauss Points
for  ii = 1:ngp
    % Define location and weight of integration
    H = H_vect(ii);
    xi = xi_vect(:, ii);

    % Find Be matrix
    [Be, det_F_Isop] = Be_6node_func(xi, [ex(1);ey(1)], [ex(2);ey(2)], [ex(3);ey(3)], [ex(4);ey(4)], [ex(5);ey(5)], [ex(6);ey(6)]);
    
    % Compute strain increment
    delta_strain = Be*(d_ae);
    % strain = strain_old + delta_strain;

    % Find Stress Trial
    stress_trial = stress_old(:, ii) + D_Mat*[delta_strain(1:2);0;delta_strain(3)];

    % Compute updated stress and updated state variables from trial stress
    [stress(:, ii), delta_strain, state(:, ii)] = mises(ptype, Mat_Prop, stress_trial', state_old(:, ii)');

    % Compute tangent material
    dstress_dstrain = dmises(ptype, Mat_Prop, stress(:, ii)', state(:, ii)');

    % Stiffness Matrix
    K_el = K_el + Be'*dstress_dstrain([1 2 4], [1 2 4])*Be*t*det_F_Isop*H;

    % Force Vector
    f_el_int = f_el_int + Be'*stress([1 2 4], ii)*t*det_F_Isop*H;

end
end