xi = sym('xi',[2,1],'real');
N1 = 0.25*(xi(1) - 1)*(xi(2) - 1); N2 = -0.25*(xi(1) + 1)*(xi(2) - 1);
N3 = 0.25*(xi(1) +1)*(xi(2) + 1); N4 = -0.25*(xi(1) - 1)*(xi(2) + 1);
%define N-matrix for the element
N_w = [N1 N2 N3 N4];
%differentiate shape functions wrt isoparam. coordinates
dN1_dxi=gradient(N1,xi);
dN2_dxi=gradient(N2,xi);
dN3_dxi=gradient(N3,xi);
dN4_dxi=gradient(N4,xi);
%introduce node positions
xe1 = sym('xe1',[2,1],'real');
xe2 = sym('xe2',[2,1],'real');
xe3 = sym('xe3',[2,1],'real');
xe4 = sym('xe4',[2,1],'real');
%introduce spatial coordinate as fcn of isoparam. coord.
x=N1*xe1+N2*xe2+N3*xe3+N4*xe4;
%compute Jacobian
Fisop=jacobian(x,xi);
detFisop = det(Fisop);
%use chain rule to compute spatial derivatives
dN1_dx=simplify( inv(Fisop)'*dN1_dxi );
dN2_dx=simplify( inv(Fisop)'*dN2_dxi );
dN3_dx=simplify( inv(Fisop)'*dN3_dxi );
dN4_dx=simplify( inv(Fisop)'*dN4_dxi );
%define B-matrix of element
Be=[dN1_dx(1), 0, dN2_dx(1), 0, dN3_dx(1), 0, dN4_dx(1), 0;
0, dN1_dx(2), 0, dN2_dx(2), 0, dN3_dx(2), 0, dN4_dx(2);
dN1_dx(2),dN1_dx(1),dN2_dx(2),dN2_dx(1),dN3_dx(2),dN3_dx(1),...
dN4_dx(2),dN4_dx(1)];
matlabFunction(Be,detFisop,N_w, 'File','Vert_Be_quad_func','Vars',{xi,xe1,xe2,xe3,xe4});