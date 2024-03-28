%% Code to Generate the functions to compute the N and B matrices for Mindlin
clc
clear 
close all

%% Mindlin plate with isoparametrical quadrilateral element 
xi = sym('xi',[2,1],'real');

% Define base functions for isoparam quadrilateral element
N1 = 0.25*(xi(1) - 1)*(xi(2) - 1); N2 = -0.25*(xi(1) + 1)*(xi(2) - 1);
N3 = 0.25*(xi(1) + 1)*(xi(2) + 1); N4 = -0.25*(xi(1) - 1)*(xi(2) + 1);

% Introduce node positions
xe1 = sym('xe1',[2,1],'real');
xe2 = sym('xe2',[2,1],'real');
xe3 = sym('xe3',[2,1],'real');
xe4 = sym('xe4',[2,1],'real');

% Introduce spatial coordinate as fcn of isoparam. coord.
x=N1*xe1+N2*xe2+N3*xe3+N4*xe4;

% Compute Jacobian
F_Isop = jacobian(x,xi);
invFisop = simplify(inv(F_Isop));
detFisop = simplify(det(F_Isop));


%% Shape Function Matrices
Ne_u = [N1 0 N2 0 N3 0 N4 0;
        0 N1 0 N2 0 N3 0 N4];

Ne_theta = Ne_u;

Ne_w =[N1 N2 N3 N4];

matlabFunction(Ne_u,Ne_theta,Ne_w, 'File','N_Mindlin_func','Vars',{xi,xe1,xe2,xe3,xe4});


%% B Matrices

% Differentiate shape functions 
dN1_dxi=gradient(N1,xi);
dN2_dxi=gradient(N2,xi);
dN3_dxi=gradient(N3,xi);
dN4_dxi=gradient(N4,xi);

% Spatial Derivatives (Chain Rule)
dN1_dx = simplify(inv(F_Isop)'*dN1_dxi );
dN2_dx = simplify(inv(F_Isop)'*dN2_dxi );
dN3_dx = simplify(inv(F_Isop)'*dN3_dxi );
dN4_dx = simplify(inv(F_Isop)'*dN4_dxi );

% B-Matrix of the element
Be_u = [dN1_dx(1),     0,     dN2_dx(1),     0,     dN3_dx(1),     0,     dN4_dx(1),     0,   ;
            0,     dN1_dx(2),     0,     dN2_dx(2),     0,     dN3_dx(2),     0,     dN4_dx(2);                                            
        dN1_dx(2), dN1_dx(1), dN2_dx(2), dN2_dx(1), dN3_dx(2), dN3_dx(1), dN4_dx(2), dN4_dx(1)];

Be_theta = Be_u;

Be_w = [dN1_dx(1), dN2_dx(1), dN3_dx(1), dN4_dx(1);
        dN1_dx(2), dN2_dx(2), dN3_dx(2), dN4_dx(2)];


matlabFunction(Be_u, Be_theta, Be_w, 'File','Be_Mindlin_func','Vars',{xi,xe1,xe2,xe3,xe4});