%% function [mesh, coord, Edof]=rectMesh(xmin, xmax, ymin, ymax, nelx, nely, nodedof)
%-----------------------------------------------------------------------------------
% INPUT:
%          xmin,xmax: min/max x-coordinates of planar plate in xy-plane
%          ymin,ymax: min/max y-coordinates of planar plate in xy-plane
%
%          nelx: number of elements along the x-direction 
%          nely: number of elements along the y-direction 
%
%          nodedof: number of dofs per node (2, 3 and 5 implemented)
%
% OUTPUT:          
%          
%          mesh: array where each column corresponds to an element and  
%                where the rows indicates the nodes associated with that 
%                element.
%
%          coord: Nnode x 2 array containing the x and y coordinates of
%                 each node (Nnode = total number of nodes)
%
%          Edof: Element topology matrix according to CALFEM
%                Each row correpsonds to an element and for that row
%                columns 2:end gives the associated degrees-of-freedom 
%                (2, 3 or 5 dofs per node)
%
%          For nodedofs = 2
%          Edof = [   1, ux1, uy1, ux2, uy2, ux3, uy3, ux4, uy4
%                     .,   .,   .,   .,   .,   .,   .,   .,  ., 
%                   nel, ux1, uy1, ux2, uy2, ux3, uy3, ux4, uy4];
%
%          For nodedofs = 3
%          Edof = [  1, w1, thetax1, thetay1, ... w4, thetax4, thetay4
%                   .,  .,       .,      ., ...  .,      .,      .
%                  nel, w1, thetax1, thetay1, ... w4, thetax4, thetay4];
%          For nodedofs = 5
%          Edof = [  1, ux1, uy1, w1, thetax1, thetay1, ... ux4, uy4, w4, thetax4, thetay4
%                   .,   .,    .,  .,      .,       .,  ...   .,   .,  .,       .,       .
%                  nel, ux1, uy1, w1, thetax1, thetay1, ... ux4, uy4, w4, thetax4, thetay4];
%----------------------------------------------------------------------------

function [mesh, coord, Edof]=rectMesh(xmin, xmax, ymin, ymax, nelx, nely, nodedof)
j=1;
k=1;
for i=1:nelx*nely
    mesh(:,i)=[(j-1)+nelx*j-(nelx-k);(j-1)+nelx*j+1-(nelx-k);(nelx+1)*(j+1)-(nelx-k);(nelx+1)*(j+1)-1-(nelx-k)];
    k=k+1;
    if (k>nelx)
        k=1;
        j=j+1;
        if(j>nely)
            break
        end
    end
end

[X,Y]=meshgrid(linspace(xmin,xmax,nelx+1),linspace(ymin,ymax,nely+1));
X=X';
Y=Y';
coord=[X(:),Y(:)];
Edof=zeros(nelx*nely,nodedof*4+1);
Edof=sparse(Edof);
Edof(:,1)=[[1:1:(nelx)*(nely)]'];
for i=1:nodedof
Edof(:,[i+1:nodedof:4*(nodedof)+1])=mesh'*nodedof-(nodedof-i);
end
