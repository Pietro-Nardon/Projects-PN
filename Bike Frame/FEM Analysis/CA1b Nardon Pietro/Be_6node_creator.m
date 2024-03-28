% Code to compute a matlab function that computes the Be matrix in a 6noded
% isotropic trangular element
 


% Define Variable
xi = sym('xi', [2, 1], 'real');

% Define shape function
N1 = (1 - xi(1) - xi(2))*(1 - 2*xi(1) - 2*xi(2));
N2 = xi(1)*(2*xi(1) - 1);
N3 = xi(2)*(2*xi(2) - 1);
N4 = 4*xi(1)*(1 - xi(1) - xi(2));
N5 = 4*xi(1)*xi(2);
N6 = 4*xi(2)*(1 - xi(1) - xi(2));

%define N-matrix for the element
Ne = [N1 0 N2 0 N3 0 N4 0 N5 0 N6 0;
      0 N1 0 N2 0 N3 0 N4 0 N5 0 N6];

% Differentiate shape functions 
dN1_dxi=gradient(N1,xi);
dN2_dxi=gradient(N2,xi);
dN3_dxi=gradient(N3,xi);
dN4_dxi=gradient(N4,xi);
dN5_dxi=gradient(N5,xi);
dN6_dxi=gradient(N6,xi);

% Introduce node position
xe1 = sym('xe1',[2,1],'real');
xe2 = sym('xe2',[2,1],'real');
xe3 = sym('xe3',[2,1],'real');
xe4 = sym('xe4',[2,1],'real');
xe5 = sym('xe5',[2,1],'real');
xe6 = sym('xe6',[2,1],'real');

% Introduce Spatial coordinates as function
x = N1*xe1 + N2*xe2 + N3*xe3 + N4*xe4 + N5*xe5 + N6*xe6;

% Compute Jacobian
F_Isop = jacobian(x, xi);
det_F_Isop = det(F_Isop);

% Spatial Derivatives (Chain Rule)
dN1_dx = simplify(inv(F_Isop)'*dN1_dxi );
dN2_dx = simplify(inv(F_Isop)'*dN2_dxi );
dN3_dx = simplify(inv(F_Isop)'*dN3_dxi );
dN4_dx = simplify(inv(F_Isop)'*dN4_dxi );
dN5_dx = simplify(inv(F_Isop)'*dN5_dxi );
dN6_dx = simplify(inv(F_Isop)'*dN6_dxi );

% B-Matrix of the element
Be = [dN1_dx(1),     0,     dN2_dx(1),     0,     dN3_dx(1),     0,     dN4_dx(1),     0,     dN5_dx(1),     0,     dN6_dx(1),     0;
          0,     dN1_dx(2),     0,     dN2_dx(2),     0,     dN3_dx(2),     0,     dN4_dx(2),     0,     dN5_dx(2),     0,     dN6_dx(2);                                            
      dN1_dx(2), dN1_dx(1), dN2_dx(2), dN2_dx(1), dN3_dx(2), dN3_dx(1), dN4_dx(2), dN4_dx(1), dN5_dx(2), dN5_dx(1), dN6_dx(2), dN6_dx(1)];


% Creation of the function
matlabFunction(Be, det_F_Isop, Ne, 'File','Be_6node_func','Vars',{xi,xe1,xe2,xe3,xe4,xe5,xe6});

