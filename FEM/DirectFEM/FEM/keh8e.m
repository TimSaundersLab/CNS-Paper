function [Ke,qe] = keh8e(Ee,Ge,Se,Xe,Set)
%***************************************************
% Keh8E:
%   Creates the element stiffness matrix of viscoelastic
%   hexahedra 8-node element.
%   Allows elastic and viscoelastic behaviour with Kelvin-Voigt (parallel
%   elastic-dahspot) and Maxwell (in series).
% Input:
%   Ee   : Strains at previous time-step
%   Ge   : element material data: [E , nu , (type)]
%          optional: type = 1 : plane stress (default)
%                    type = 2 : plane strain
%   Se   : Stresses at previous time-step.Kelvin: elastic stresses
%   Xe   : coordinates Xe = [x1 y1; x2 y2; x3 y3; x4 y4]
%   Set.Rheology='Elastic', 'Kelvin', 'Maxwell'. Viscosity as a function of full strains (deviatoric and non-deviatoric).
% Output:
%   Ke   : element stiffness matrix.
%   qe   : force vector. -qe is sent to rhs.
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
if strcmp(Set.Rheology(1:6),'Maxwel')==1
    DMaxwell=inv(inv(D)+Set.dt/Set.eta*Set.theta*eye(6)); % =D*inv(eye(6)+Set.dt/Set.eta*Set.theta)
    qMaxwell=inv(D)-Set.dt/Set.eta*(1-Set.theta)*eye(6);
end
   
% determine number of nodes per element
nnodes = size(Xe,1);

% Initialize stiffness matrix.
Ke = zeros(3*nnodes);
qe=zeros(3*nnodes,1);
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
            
            if strcmp(Set.Rheology(1:6),'Kelvin')==1
                Sgp(1:6)=Se(1,gp,:);
                Egp(1:6)=Ee(1,gp,:);
                Ke = Ke + Set.theta*w(i)*w(j)*w(k)*( B'*D*B )*det(Jt)...
                +Set.eta/Set.dt*w(i)*w(j)*w(k)*( B'*B)*det(Jt);
                qe = qe + (1-Set.theta)*w(i)*w(j)*w(k)*( B'*Sgp')*det(Jt)...
                    -Set.eta/Set.dt*w(i)*w(j)*w(k)*( B'*Egp')*det(Jt);
            elseif strcmp(Set.Rheology(1:6),'Maxwel')==1
                Sgp(1:6)=Se(1,gp,:);
                Egp(1:6)=Ee(1,gp,:);
                Ke = Ke + w(i)*w(j)*w(k)*( B'*DMaxwell*B )*det(Jt);
                qe = qe + w(i)*w(j)*w(k)*...
                    B'*DMaxwell*(qMaxwell*Sgp'-Egp')*det(Jt);
            else
                Ke = Ke + w(i)*w(j)*w(k)*( B'*D*B )*det(Jt);
            end
            
        end
    end
    
end