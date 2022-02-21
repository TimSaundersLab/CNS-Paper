function [C,X,Bij]=GenerateMesh(dx,dy,dz,nx,ny,nz)
% USAGE:
%  [X,C]=GenerateMesh(dx,dy,dz,nx,ny,nz)
% INPUT:
% nx  = number of elements (divisions) along x (Horizontal)
% dx = size of each element along x
% ny = number of elements (divisions) along y (Vertical)
% dx = size of each element along y
% nz = number of elements (divisions) along z
% dz = size of each element along z (if dz==0, 2D)
% OUTPUT:
% X(:,i)    =nodal coordiantes of node i
% C(:,e) =nodal connectivities of element e
% Bij    = matrix with position of nodes at each grid point
%      2D: Bij(1,1) node at bottom left,
%          Bij(1,ny+1) node at top left,
%          Bij(nx+1,ny+1) node at top right
%
%           ^            +----> ny
%           | coord(X)   |
%        ny |            |nx   Bij
%           +---> nx     v
%      3D: Bij(1,1,1)          node at bottom left on z=0,
%          Bij(1,ny+1,1)       node at top left on z=0,
%          Bij(nx+1,ny+1,nz+1) node at top right on z=nz*dx
if ~exist('nz','var')
    nz=0;
end
if isempty(dz)
    dz=0;
end
if abs(dy)<eps || ny==0 % 1D
    X=zeros(nx+1,1);
    C=zeros(nx,2);
    for i=1:nx+1 % Horizontal
        X(i)=(i-1)*dx;
        if i~=nx+1
            C(i,:)=[i i+1];
        end
    end
    Bij=1:nx+1;
elseif abs(dz)<eps || nz<eps % 2D
    if nargout >2
        Bij=zeros(nx+1,ny+1);
    end
    X=zeros((nx+1)*(ny+1),2);
    C=zeros(nx*ny,4);
    for j=1:ny+1 % horizontal (left to right)
        for i=1:nx+1 % Vertical (bottom to top)
            nn=(i-1)*(ny+1)+j;
            X(nn,:)=[(i-1)*dx,(j-1)*dy];
            if nargout>2
                Bij(i,j)=nn;
            end
            if i~=nx+1 && j~=ny+1
                ne=(nx-i)*ny+j;
                nbl=nn;
                bot=[nbl,nbl+1,nbl+ny+2,nbl+ny+1];
                C(ne,:)=bot;
            end
        end
    end
else
    if nargout>2
        Bij=zeros(nx+1,ny+1,nz+1);
    end
    X=zeros((nx+1)*(ny+1)*(nz+1),3);
    C=zeros(nx*ny*nz,8);
    for k=1:nz+1
        for j=1:ny+1 % horizontal (left to right)
            for i=1:nx+1 % Vertical (bottom to top)
                nn=(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i;
                X(nn,:)=[(i-1)*dx,(j-1)*dy,(k-1)*dz];
                if nargout>2
                    Bij(i,j,k)=nn;
                end
                if k~=nz+1 && i~=nx+1 && j~=ny+1
                    ne=(k-1)*nx*ny+(j-1)*nx+i;
                    nbl=nn;
                    bot=[nbl,nbl+1,nbl+nx+2,nbl+nx+1];
                    top=bot+(ny+1)*(nx+1);
                    C(ne,:)=[bot,top];
                end
            end
        end
    end
end
% Add material number
C=[C ones(size(C,1),1)];
