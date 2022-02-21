function [qe,Ee,Se,Sev] = qeh8e(En,Ge,Ue,Sn,Set,Xe)
%***************************************************
% Keh8E:
%   Creates the element elastic forcee vector for
%   hexahedra 8-node element.
% Syntax:
%   Ke = qeh8e(Xe,Ge)
% Input:
%   Xe   : coordinates Xe = [x1 y1; x2 y2; x3 y3; x4 y4]
%   Ge   : element material data: [E , nu , (type)]
%   Een  : previous time-step elemental strain
%   Sen  : previous time-step elemental stress
%   Ue   : element displacement vector: 
%               Ue = [u1;v1;w1;u2;v2;w2;...u8;v8;w8]
% Output:
%   Ke   : element stiffness matrix.
%   Ee   : current elemental stress
%   Se   : current elastic elemental strains
%   Sev  : Kelvin: current viscous stresses
%          Maxwell: 0 (all stresses included in elastic part)
% Date:
%   Version 1.0    04.05.2014
%***************************************************

% Gauss abscissae and weights.
r = [-1 1]/sqrt(3);
w = [ 1 1];

% Set isotropic elasticity matrix
E  = Ge(1);
nu = Ge(2);
f=E/(1-2*nu)/(1+nu);

D(1:3,1:3)  = f* ...
    [ 1-nu  nu     nu
      nu    1-nu   nu
      nu    nu     1-nu];
D(4,4)=f*(1-2*nu)/2;
D(5,5)=D(4,4);
D(6,6)=D(5,5);

% determine number of nodes per element
nnodes = size(Xe,1);
Ue=reshape(Ue',[3*nnodes,1]);

if strcmp(Set.Rheology(1:6),'Maxwel')==1
    iMaxwell=inv(eye(6)+Set.dt/Set.eta*Set.theta*D);
    qMaxwell=eye(6)-Set.dt/Set.eta*(1-Set.theta)*D;
end
% Initialize internal force vector.
qe = zeros(3*nnodes,1);
Ee=zeros(8,6);
Se=zeros(8,6);
Sev=zeros(8,6);
zerosS=zeros(size(Se));
% Gauss integration of stiffness matrix.
for i = 1:2
    for j = 1:2
        for k = 1:2
            % organize the Gauss points as the element nodes
            gp=i+3*(j-1) - 2*(i-1)*(j-1)+4*(k-1);
            % Parametric derivatives:
            dN(1:3,1:4) = [ -(1-r(j))  (1-r(j))  (1+r(j)) -(1+r(j))
                            -(1-r(i)) -(1+r(i))  (1+r(i))  (1-r(i))
                (1-r(i))*(1-r(j)) (1+r(i))*(1-r(j))  (1+r(i))*(1+r(j)) (1-r(i))*(1+r(j))]/8;
            dN(3,5:8) = dN(3,1:4);
            dN(3,1:4) = -dN(3,1:4);
            dN(1:2,5:8) = dN(1:2,1:4)*(1+r(k));
            dN(1:2,1:4) = dN(1:2,1:4)*(1-r(k));
            % transform to global coordinates
            Jt = dN*Xe;
            dN = Jt\dN;
            
            % set up 4 node part of the gradient matrix
            B=zeros(6,3*nnodes);
            B(1,1:3:24)  =dN(1,:);
            B(2,2:3:24)  =dN(2,:);
            B(3,3:3:24)  =dN(3,:);
            B(4,1:3:24)  =dN(2,:);
            B(4,2:3:24)  =dN(1,:);
            B(5,1:3:24)  =dN(3,:);
            B(5,3:3:24)  =dN(1,:);
            B(6,2:3:24)  =dN(3,:);
            B(6,3:3:24)  =dN(2,:);
            
            Ee(gp,:)=(B*Ue)'; % Total Strains in Maxwell
            if strcmp(Set.Rheology(1:6),'Kelvin')==1
                Se(gp,:)=Ee(gp,:)*D; % Elastic stresses in Kelvin viscoelasticity
                Sev(gp,:)=Set.eta/Set.dt*(Ee(gp,:)-En(gp,:));% Viscous stresses
                qe = qe + w(i)*w(j)*w(k)*( B'*(Se(gp,:)+Sev(gp,:))')*det(Jt);
            elseif strcmp(Set.Rheology(1:6),'Maxwel')==1
                Se(gp,:)=iMaxwell*(D*(Ee(gp,:)'-En(gp,:)')+qMaxwell*Sn(gp,:)'); 
                Sev(gp,:)=Se(gp,:);
                qe = qe + w(i)*w(j)*w(k)* B'*Se(gp,:)'*det(Jt);
            else % Elastic
                Se(gp,:)=Ee(gp,:)*D; 
                Sev(gp,:)=zeros(1,size(Se,2));
                qe = qe + w(i)*w(j)*w(k)*( B'*(Se(gp,:)+Sev(gp,:))')*det(Jt);
            end
            
        end
    end
    
end